#!/usr/bin/env python3
"""Apply manual MedDRA overrides for unresolved reaction terms in Step5 outputs."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import sqlite3
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Apply manual reaction MedDRA overrides")
    p.add_argument("--step2-db", type=Path, default=Path("outputs/step2/faers_step2.sqlite"))
    p.add_argument("--step5-db", type=Path, default=Path("outputs/step5/faers_step5.sqlite"))
    p.add_argument(
        "--meddra-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "reference_data" / "MedDRA" / "MedDRA_29_0_English",
        help="Path containing MedAscii/pt.asc for fallback PT code->name lookup",
    )
    p.add_argument(
        "--overrides-csv",
        type=Path,
        default=Path("resources/step5_manual/reaction_term_manual_overrides_2026-04-13.csv"),
    )
    p.add_argument("--output-dir", type=Path, default=Path("outputs/step5"))
    return p.parse_args()


def scalar(conn: sqlite3.Connection, sql: str) -> int:
    return int(conn.execute(sql).fetchone()[0])


def load_overrides(path: Path) -> list[tuple[str, str, str]]:
    rows: list[tuple[str, str, str]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for r in reader:
            term = (r.get("pt_norm") or "").strip().upper()
            code = (r.get("pt_meddra_code") or "").strip()
            reason = (r.get("override_reason") or "").strip()
            if term and code:
                rows.append((term, code, reason))
    return rows


def find_medascii_dir(meddra_dir: Path) -> Path | None:
    if not meddra_dir.exists():
        return None
    direct = meddra_dir / "MedAscii"
    if direct.exists() and direct.is_dir():
        return direct
    for d in meddra_dir.rglob("MedAscii"):
        if d.is_dir():
            return d
    return None


def load_pt_code_to_name_from_meddra(meddra_dir: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    medascii = find_medascii_dir(meddra_dir)
    if medascii is None:
        return out
    pt_path = medascii / "pt.asc"
    if not pt_path.exists():
        return out
    with pt_path.open("r", encoding="latin-1", errors="replace", newline="") as f:
        for line in f:
            row = line.rstrip("\r\n").split("$")
            if len(row) < 2:
                continue
            code = row[0].strip()
            name = row[1].strip()
            if code and name:
                out[code] = name
    return out


def rebuild_step5_report(conn: sqlite3.Connection, output_dir: Path, started_at: str, finished_at: str) -> None:
    mapping_source_counts = {
        r[0]: int(r[1])
        for r in conn.execute(
            """
            SELECT pt_mapping_source, COUNT(*)
            FROM reaction_term_meddra_map
            GROUP BY pt_mapping_source
            ORDER BY COUNT(*) DESC, pt_mapping_source
            """
        )
    }

    total_rows = scalar(conn, "SELECT COUNT(*) FROM s2.case_reaction")
    mapped_rows = scalar(
        conn,
        """
        SELECT COUNT(*)
        FROM s2.case_reaction c
        INNER JOIN reaction_term_meddra_map m ON c.pt_norm = m.pt_norm
        WHERE TRIM(COALESCE(m.pt_meddra_code,'')) <> ''
        """,
    )
    total_terms = scalar(conn, "SELECT COUNT(*) FROM reaction_term_meddra_map")
    mapped_terms = scalar(
        conn,
        "SELECT COUNT(*) FROM reaction_term_meddra_map WHERE TRIM(COALESCE(pt_meddra_code,'')) <> ''",
    )

    report = {
        "started_at": started_at,
        "finished_at": finished_at,
        "meddra_release": "29.0$English$$$$",
        "step2_db": "outputs\\step2\\faers_step2.sqlite",
        "step5_db": "outputs\\step5\\faers_step5.sqlite",
        "reaction_rows_total": total_rows,
        "reaction_rows_mapped": mapped_rows,
        "reaction_rows_unmapped": total_rows - mapped_rows,
        "reaction_row_coverage_mapped": (mapped_rows / total_rows) if total_rows else 0.0,
        "term_total": total_terms,
        "term_mapped": mapped_terms,
        "term_unmapped": total_terms - mapped_terms,
        "term_coverage_mapped": (mapped_terms / total_terms) if total_terms else 0.0,
        "pt_distinct_mapped": scalar(conn, "SELECT COUNT(DISTINCT pt_code) FROM reaction_base_meddra"),
        "primaryid_distinct_mapped": scalar(conn, "SELECT COUNT(DISTINCT primaryid) FROM reaction_base_meddra"),
        "mapping_source_counts": mapping_source_counts,
        "exports": {
            "reaction_term_meddra_map_rows": scalar(conn, "SELECT COUNT(*) FROM reaction_term_meddra_map"),
            "reaction_unresolved_rows": scalar(conn, "SELECT COUNT(*) FROM reaction_unresolved_top"),
            "reaction_unresolved_top1000_rows": min(1000, scalar(conn, "SELECT COUNT(*) FROM reaction_unresolved_top")),
            "pt_case_counts_rows": scalar(conn, "SELECT COUNT(*) FROM pt_case_counts_meddra"),
        },
    }
    (output_dir / "step5_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()
    started = dt.datetime.now().isoformat(timespec="seconds")

    if not args.step2_db.exists():
        print(f"[ERROR] step2 db not found: {args.step2_db}", flush=True)
        return 1
    if not args.step5_db.exists():
        print(f"[ERROR] step5 db not found: {args.step5_db}", flush=True)
        return 1
    if not args.overrides_csv.exists():
        print(f"[ERROR] overrides csv not found: {args.overrides_csv}", flush=True)
        return 1

    overrides = load_overrides(args.overrides_csv)
    if not overrides:
        print("[ERROR] no valid overrides loaded", flush=True)
        return 1

    fallback_code_name = load_pt_code_to_name_from_meddra(args.meddra_dir)

    conn = sqlite3.connect(args.step5_db)
    try:
        conn.execute("PRAGMA journal_mode = WAL;")
        conn.execute("PRAGMA synchronous = NORMAL;")
        conn.execute("PRAGMA temp_store = FILE;")
        conn.execute("ATTACH DATABASE ? AS s2", (str(args.step2_db),))

        before_unmapped_terms = scalar(
            conn,
            "SELECT COUNT(*) FROM reaction_term_meddra_map WHERE TRIM(COALESCE(pt_meddra_code,'')) = ''",
        )
        before_mapped_rows = scalar(
            conn,
            """
            SELECT COUNT(*)
            FROM s2.case_reaction c
            INNER JOIN reaction_term_meddra_map m ON c.pt_norm = m.pt_norm
            WHERE TRIM(COALESCE(m.pt_meddra_code,'')) <> ''
            """,
        )

        conn.executescript(
            """
            CREATE TABLE IF NOT EXISTS reaction_term_manual_overrides (
                pt_norm TEXT PRIMARY KEY,
                pt_meddra_code TEXT NOT NULL,
                override_reason TEXT,
                applied_at TEXT
            );
            DELETE FROM reaction_term_manual_overrides;
            """
        )

        applied = 0
        missing_terms = []
        unresolved_code_names = []
        applied_at = dt.datetime.now().isoformat(timespec="seconds")

        for term, code, reason in overrides:
            exists = conn.execute(
                "SELECT 1 FROM reaction_term_meddra_map WHERE pt_norm = ? LIMIT 1",
                (term,),
            ).fetchone()
            if not exists:
                missing_terms.append(term)
                continue

            name_row = conn.execute(
                """
                SELECT pt_meddra_pt
                FROM reaction_term_meddra_map
                WHERE pt_meddra_code = ? AND TRIM(COALESCE(pt_meddra_pt,'')) <> ''
                LIMIT 1
                """,
                (code,),
            ).fetchone()
            if not name_row:
                pt_name = fallback_code_name.get(code, "")
                if not pt_name:
                    unresolved_code_names.append((term, code))
                    continue
            else:
                pt_name = name_row[0]
            conn.execute(
                """
                UPDATE reaction_term_meddra_map
                SET pt_meddra_code = ?,
                    pt_meddra_pt = ?,
                    pt_mapping_source = 'MANUAL_TERM_OVERRIDE'
                WHERE pt_norm = ?
                """,
                (code, pt_name, term),
            )
            conn.execute(
                """
                INSERT INTO reaction_term_manual_overrides(pt_norm, pt_meddra_code, override_reason, applied_at)
                VALUES (?, ?, ?, ?)
                ON CONFLICT(pt_norm) DO UPDATE SET
                    pt_meddra_code=excluded.pt_meddra_code,
                    override_reason=excluded.override_reason,
                    applied_at=excluded.applied_at
                """,
                (term, code, reason, applied_at),
            )
            applied += 1

        conn.commit()

        override_terms = [t for t, _, _ in overrides]
        q_marks = ",".join("?" for _ in override_terms)
        if override_terms:
            conn.execute(
                f"""
                INSERT INTO reaction_base_meddra(primaryid, caseid, pt_code, pt_name, map_source)
                SELECT DISTINCT
                    c.primaryid,
                    c.caseid,
                    m.pt_meddra_code AS pt_code,
                    m.pt_meddra_pt AS pt_name,
                    m.pt_mapping_source AS map_source
                FROM s2.case_reaction c
                INNER JOIN reaction_term_meddra_map m ON c.pt_norm = m.pt_norm
                WHERE c.pt_norm IN ({q_marks})
                  AND TRIM(COALESCE(m.pt_meddra_code,'')) <> ''
                """,
                override_terms,
            )
            conn.commit()

        conn.executescript(
            """
            DROP TABLE IF EXISTS reaction_base_meddra_dedup;
            CREATE TABLE reaction_base_meddra_dedup AS
            SELECT DISTINCT primaryid, caseid, pt_code, pt_name, map_source
            FROM reaction_base_meddra;
            DROP TABLE reaction_base_meddra;
            ALTER TABLE reaction_base_meddra_dedup RENAME TO reaction_base_meddra;
            CREATE INDEX idx_reaction_base_meddra_pid ON reaction_base_meddra(primaryid);
            CREATE INDEX idx_reaction_base_meddra_ptcode ON reaction_base_meddra(pt_code);
            CREATE INDEX idx_reaction_base_meddra_ptname ON reaction_base_meddra(pt_name);

            DROP TABLE IF EXISTS reaction_unresolved_top;
            CREATE TABLE reaction_unresolved_top AS
            SELECT
                c.pt_norm,
                COUNT(*) AS cnt
            FROM s2.case_reaction c
            INNER JOIN reaction_term_meddra_map m ON c.pt_norm = m.pt_norm
            WHERE TRIM(COALESCE(m.pt_meddra_code,'')) = ''
            GROUP BY c.pt_norm
            ORDER BY cnt DESC, c.pt_norm;
            CREATE INDEX idx_unresolved_cnt ON reaction_unresolved_top(cnt);

            DROP TABLE IF EXISTS pt_case_counts_meddra;
            CREATE TABLE pt_case_counts_meddra AS
            SELECT
                pt_code,
                pt_name,
                COUNT(*) AS n_case_reaction_rows,
                COUNT(DISTINCT primaryid) AS n_primaryid
            FROM reaction_base_meddra
            GROUP BY pt_code, pt_name
            ORDER BY n_primaryid DESC, pt_code;
            CREATE INDEX idx_pt_case_counts_code ON pt_case_counts_meddra(pt_code);
            CREATE INDEX idx_pt_case_counts_npid ON pt_case_counts_meddra(n_primaryid);
            """
        )
        conn.commit()

        after_unmapped_terms = scalar(
            conn,
            "SELECT COUNT(*) FROM reaction_term_meddra_map WHERE TRIM(COALESCE(pt_meddra_code,'')) = ''",
        )
        after_mapped_rows = scalar(
            conn,
            """
            SELECT COUNT(*)
            FROM s2.case_reaction c
            INNER JOIN reaction_term_meddra_map m ON c.pt_norm = m.pt_norm
            WHERE TRIM(COALESCE(m.pt_meddra_code,'')) <> ''
            """,
        )

        args.output_dir.mkdir(parents=True, exist_ok=True)
        with (args.output_dir / "reaction_unresolved_top1000.csv").open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["pt_norm", "cnt"])
            for r in conn.execute(
                "SELECT pt_norm, cnt FROM reaction_unresolved_top ORDER BY cnt DESC, pt_norm LIMIT 1000"
            ):
                w.writerow(r)

        with (args.output_dir / "reaction_unresolved_top.csv").open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["pt_norm", "cnt"])
            for r in conn.execute("SELECT pt_norm, cnt FROM reaction_unresolved_top ORDER BY cnt DESC, pt_norm"):
                w.writerow(r)

        with (args.output_dir / "reaction_term_meddra_map.csv").open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["pt_norm", "pt_manual_code_hint", "pt_meddra_code", "pt_meddra_pt", "pt_mapping_source"])
            for r in conn.execute(
                """
                SELECT pt_norm, pt_manual_code_hint, pt_meddra_code, pt_meddra_pt, pt_mapping_source
                FROM reaction_term_meddra_map
                ORDER BY pt_mapping_source, pt_norm
                """
            ):
                w.writerow(r)

        with (args.output_dir / "pt_case_counts_meddra.csv").open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["pt_code", "pt_name", "n_case_reaction_rows", "n_primaryid"])
            for r in conn.execute(
                "SELECT pt_code, pt_name, n_case_reaction_rows, n_primaryid FROM pt_case_counts_meddra"
            ):
                w.writerow(r)

        finished = dt.datetime.now().isoformat(timespec="seconds")
        report = {
            "started_at": started,
            "finished_at": finished,
            "overrides_csv": str(args.overrides_csv),
            "override_candidates": len(overrides),
            "override_applied": applied,
            "missing_terms_in_map": missing_terms,
            "codes_without_name_lookup": unresolved_code_names,
            "before_unmapped_terms": before_unmapped_terms,
            "after_unmapped_terms": after_unmapped_terms,
            "before_mapped_rows": before_mapped_rows,
            "after_mapped_rows": after_mapped_rows,
            "rows_gained_by_override": after_mapped_rows - before_mapped_rows,
        }
        (args.output_dir / "step5_manual_override_report.json").write_text(
            json.dumps(report, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

        rebuild_step5_report(conn, args.output_dir, started, finished)
        conn.execute("PRAGMA wal_checkpoint(TRUNCATE);")

        print("[DONE] manual overrides applied", flush=True)
        print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    raise SystemExit(main())
