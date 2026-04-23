#!/usr/bin/env python3
"""Step 4 supplementary export from existing InChIKey tables.

Purpose:
- continue Step4 from previously materialized structure tables
- rebuild PT-level aggregate exports without repeating upstream normalization
- provide a recovery path after interrupted long-running aggregations

Inputs:
- outputs/step4/faers_step4.sqlite

Outputs:
- outputs/step4/*.csv
- outputs/step4/*.json
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import sqlite3
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Step 4: finalize PT-level exports from existing InChIKey tables")
    p.add_argument("--output-dir", type=Path, default=Path("outputs/step4"))
    p.add_argument("--min-pair-cases", type=int, default=20)
    p.add_argument("--include-any", action="store_true")
    return p.parse_args()


def scalar(conn: sqlite3.Connection, q: str) -> int:
    return int(conn.execute(q).fetchone()[0])


def table_exists(conn: sqlite3.Connection, table_name: str) -> bool:
    q = "SELECT 1 FROM sqlite_master WHERE type='table' AND name=? LIMIT 1"
    return conn.execute(q, (table_name,)).fetchone() is not None


def export_csv(conn: sqlite3.Connection, sql: str, out_path: Path) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cur = conn.cursor()
    cur.execute(sql)
    cols = [d[0] for d in cur.description]
    n = 0
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for row in cur:
            w.writerow(row)
            n += 1
    return n


def main() -> int:
    args = parse_args()
    started = dt.datetime.now().isoformat(timespec="seconds")
    step4_db = args.output_dir / "faers_step4.sqlite"
    if not step4_db.exists():
        print(f"[ERROR] step4 db not found: {step4_db}", flush=True)
        return 1

    conn = sqlite3.connect(step4_db)
    try:
        conn.execute("PRAGMA journal_mode = WAL;")
        conn.execute("PRAGMA synchronous = NORMAL;")
        conn.execute("PRAGMA temp_store = FILE;")

        required = [
            "term_inchi_raw",
            "inchi_norm_map",
            "term_inchikey_map",
            "inchikey_meta",
            "reaction_base",
            "case_inchikey_any",
            "case_inchikey_psss",
            "inchikey_case_counts_any",
            "inchikey_case_counts_psss",
        ]
        missing = [t for t in required if not table_exists(conn, t)]
        if missing:
            print(f"[ERROR] missing prerequisite tables: {missing}", flush=True)
            return 1

        if args.include_any:
            print("[INFO] build inchikey_pt_counts_any ...", flush=True)
            conn.executescript(
                """
                DROP TABLE IF EXISTS inchikey_pt_counts_any;
                CREATE TABLE inchikey_pt_counts_any AS
                SELECT
                    c.inchikey,
                    r.pt_norm,
                    COUNT(*) AS n_cases
                FROM case_inchikey_any c
                INNER JOIN reaction_base r ON r.primaryid = c.primaryid
                GROUP BY c.inchikey, r.pt_norm;
                CREATE INDEX idx_ik_pt_any_key ON inchikey_pt_counts_any(inchikey);
                CREATE INDEX idx_ik_pt_any_pt ON inchikey_pt_counts_any(pt_norm);
                """
            )
            conn.commit()
        else:
            conn.executescript(
                """
                DROP TABLE IF EXISTS inchikey_pt_counts_any;
                CREATE TABLE inchikey_pt_counts_any (
                    inchikey TEXT,
                    pt_norm TEXT,
                    n_cases INTEGER
                );
                """
            )
            conn.commit()

        print("[INFO] build inchikey_pt_counts_psss ...", flush=True)
        conn.executescript(
            """
            DROP TABLE IF EXISTS inchikey_pt_counts_psss;
            CREATE TABLE inchikey_pt_counts_psss AS
            SELECT
                c.inchikey,
                r.pt_norm,
                COUNT(*) AS n_cases
            FROM case_inchikey_psss c
            INNER JOIN reaction_base r ON r.primaryid = c.primaryid
            GROUP BY c.inchikey, r.pt_norm;
            CREATE INDEX idx_ik_pt_psss_key ON inchikey_pt_counts_psss(inchikey);
            CREATE INDEX idx_ik_pt_psss_pt ON inchikey_pt_counts_psss(pt_norm);
            """
        )
        conn.commit()

        print("[INFO] build model/index/matrix ...", flush=True)
        conn.execute("DROP TABLE IF EXISTS inchikey_pt_psss_model")
        conn.execute(
            """
            CREATE TABLE inchikey_pt_psss_model AS
            SELECT inchikey, pt_norm, n_cases
            FROM inchikey_pt_counts_psss
            WHERE n_cases >= ?
            """,
            (args.min_pair_cases,),
        )
        conn.commit()

        conn.executescript(
            """
            DROP TABLE IF EXISTS inchikey_index;
            CREATE TABLE inchikey_index AS
            SELECT
                ROW_NUMBER() OVER (ORDER BY inchikey) AS inchikey_id,
                inchikey
            FROM (SELECT DISTINCT inchikey FROM inchikey_pt_psss_model);
            CREATE INDEX idx_inchikey_index_key ON inchikey_index(inchikey);

            DROP TABLE IF EXISTS pt_index;
            CREATE TABLE pt_index AS
            SELECT
                ROW_NUMBER() OVER (ORDER BY pt_norm) AS pt_id,
                pt_norm
            FROM (SELECT DISTINCT pt_norm FROM inchikey_pt_psss_model);
            CREATE INDEX idx_pt_index_pt ON pt_index(pt_norm);

            DROP TABLE IF EXISTS inchikey_pt_matrix_psss;
            CREATE TABLE inchikey_pt_matrix_psss AS
            SELECT
                i.inchikey_id,
                p.pt_id,
                m.n_cases
            FROM inchikey_pt_psss_model m
            INNER JOIN inchikey_index i ON i.inchikey = m.inchikey
            INNER JOIN pt_index p ON p.pt_norm = m.pt_norm;
            CREATE INDEX idx_matrix_ik ON inchikey_pt_matrix_psss(inchikey_id);
            CREATE INDEX idx_matrix_pt ON inchikey_pt_matrix_psss(pt_id);
            """
        )
        conn.commit()

        print("[INFO] export csv/report ...", flush=True)
        n_inchi_norm = export_csv(
            conn,
            "SELECT * FROM inchi_norm_map ORDER BY status, inchi_raw",
            args.output_dir / "inchi_normalization_map.csv",
        )
        n_term_map = export_csv(
            conn,
            "SELECT * FROM term_inchikey_map ORDER BY n_rows DESC, term_std",
            args.output_dir / "term_inchikey_map.csv",
        )
        n_meta = export_csv(
            conn,
            "SELECT * FROM inchikey_meta ORDER BY term_count DESC, inchikey",
            args.output_dir / "inchikey_meta.csv",
        )
        n_case_any = export_csv(
            conn,
            "SELECT * FROM inchikey_case_counts_any ORDER BY n_cases DESC, inchikey",
            args.output_dir / "inchikey_case_counts_any.csv",
        )
        n_case_psss = export_csv(
            conn,
            "SELECT * FROM inchikey_case_counts_psss ORDER BY n_cases DESC, inchikey",
            args.output_dir / "inchikey_case_counts_psss.csv",
        )
        n_pair_psss = export_csv(
            conn,
            "SELECT * FROM inchikey_pt_psss_model ORDER BY n_cases DESC, inchikey, pt_norm",
            args.output_dir / f"inchikey_pt_counts_psss_min{args.min_pair_cases}.csv",
        )
        n_ik_idx = export_csv(
            conn,
            "SELECT * FROM inchikey_index ORDER BY inchikey_id",
            args.output_dir / "inchikey_index.csv",
        )
        n_pt_idx = export_csv(
            conn,
            "SELECT * FROM pt_index ORDER BY pt_id",
            args.output_dir / "pt_index.csv",
        )
        n_edges = export_csv(
            conn,
            "SELECT * FROM inchikey_pt_matrix_psss ORDER BY inchikey_id, pt_id",
            args.output_dir / f"inchikey_pt_matrix_psss_min{args.min_pair_cases}.csv",
        )
        export_csv(
            conn,
            "SELECT * FROM inchikey_pt_counts_any ORDER BY n_cases DESC, inchikey, pt_norm LIMIT 20000",
            args.output_dir / "inchikey_pt_counts_any_top20000.csv",
        )
        export_csv(
            conn,
            "SELECT * FROM inchikey_pt_counts_psss ORDER BY n_cases DESC, inchikey, pt_norm LIMIT 20000",
            args.output_dir / "inchikey_pt_counts_psss_top20000.csv",
        )

        finished = dt.datetime.now().isoformat(timespec="seconds")
        summary = {
            "started_at": started,
            "finished_at": finished,
            "runner": "step4_finalize_pt_export.py",
            "min_pair_cases": args.min_pair_cases,
            "include_any": bool(args.include_any),
            "term_inchi_raw_rows": scalar(conn, "SELECT COUNT(*) FROM term_inchi_raw"),
            "term_inchi_raw_distinct_terms": scalar(conn, "SELECT COUNT(DISTINCT term_std) FROM term_inchi_raw"),
            "term_inchi_raw_distinct_inchi": scalar(conn, "SELECT COUNT(DISTINCT inchi_raw) FROM term_inchi_raw"),
            "raw_inchi_total": scalar(conn, "SELECT COUNT(*) FROM inchi_norm_map"),
            "raw_inchi_norm_ok": scalar(conn, "SELECT COUNT(*) FROM inchi_norm_map WHERE status='OK'"),
            "inchi_norm_rows": n_inchi_norm,
            "term_inchikey_map_rows": n_term_map,
            "inchikey_meta_rows": n_meta,
            "reaction_base_rows": scalar(conn, "SELECT COUNT(*) FROM reaction_base"),
            "case_inchikey_any_rows": scalar(conn, "SELECT COUNT(*) FROM case_inchikey_any"),
            "case_inchikey_psss_rows": scalar(conn, "SELECT COUNT(*) FROM case_inchikey_psss"),
            "inchikey_case_any_rows": n_case_any,
            "inchikey_case_psss_rows": n_case_psss,
            "inchikey_pt_any_rows": scalar(conn, "SELECT COUNT(*) FROM inchikey_pt_counts_any"),
            "inchikey_pt_psss_rows": scalar(conn, "SELECT COUNT(*) FROM inchikey_pt_counts_psss"),
            "inchikey_pt_psss_model_rows": n_pair_psss,
            "inchikey_index_rows": n_ik_idx,
            "pt_index_rows": n_pt_idx,
            "matrix_edges_rows": n_edges,
            "max_pair_cases_any": scalar(conn, "SELECT COALESCE(MAX(n_cases),0) FROM inchikey_pt_counts_any"),
            "max_pair_cases_psss": scalar(conn, "SELECT COALESCE(MAX(n_cases),0) FROM inchikey_pt_counts_psss"),
            "inchikey_status_breakdown": {
                "ok": scalar(conn, "SELECT COUNT(*) FROM inchi_norm_map WHERE status='OK'"),
                "failed": scalar(conn, "SELECT COUNT(*) FROM inchi_norm_map WHERE status<>'OK'"),
                "multi_fragment_before": scalar(conn, "SELECT COUNT(*) FROM inchi_norm_map WHERE n_frag_before > 1"),
            },
        }
        (args.output_dir / "step4_report.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

        conn.execute("PRAGMA wal_checkpoint(TRUNCATE);")
        print("[DONE] step4 finalize completed", flush=True)
        print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    raise SystemExit(main())
