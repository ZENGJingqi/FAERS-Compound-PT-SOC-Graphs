#!/usr/bin/env python3
"""Build MedDRA-standardized InChIKey-PT tables using Step4 + Step5 outputs."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import sqlite3
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build InChIKey-PT (MedDRA standardized, PS/SS)")
    p.add_argument("--step4-db", type=Path, default=Path("outputs/step4/faers_step4.sqlite"))
    p.add_argument("--step5-db", type=Path, default=Path("outputs/step5/faers_step5.sqlite"))
    p.add_argument("--output-dir", type=Path, default=Path("outputs/step5"))
    p.add_argument("--min-pair-cases", type=int, default=20)
    p.add_argument("--rebuild", action="store_true")
    return p.parse_args()


def scalar(conn: sqlite3.Connection, sql: str) -> int:
    return int(conn.execute(sql).fetchone()[0])


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

    if not args.step4_db.exists():
        print(f"[ERROR] step4 db not found: {args.step4_db}", flush=True)
        return 1
    if not args.step5_db.exists():
        print(f"[ERROR] step5 db not found: {args.step5_db}", flush=True)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(args.step5_db)
    try:
        conn.execute("PRAGMA journal_mode = WAL;")
        conn.execute("PRAGMA synchronous = NORMAL;")
        conn.execute("PRAGMA temp_store = FILE;")
        conn.execute("ATTACH DATABASE ? AS s4", (str(args.step4_db),))

        if args.rebuild:
            conn.executescript(
                """
                DROP TABLE IF EXISTS inchikey_case_counts_psss_meddra;
                DROP TABLE IF EXISTS inchikey_pt_counts_psss_meddra;
                DROP TABLE IF EXISTS inchikey_pt_psss_meddra_model;
                DROP TABLE IF EXISTS inchikey_index_meddra;
                DROP TABLE IF EXISTS pt_index_meddra;
                DROP TABLE IF EXISTS inchikey_pt_matrix_psss_meddra;
                """
            )
            conn.commit()

        print("[INFO] building inchikey_case_counts_psss_meddra ...", flush=True)
        conn.executescript(
            """
            DROP TABLE IF EXISTS inchikey_case_counts_psss_meddra;
            CREATE TABLE inchikey_case_counts_psss_meddra AS
            SELECT inchikey, COUNT(*) AS n_cases
            FROM s4.case_inchikey_psss
            GROUP BY inchikey;
            CREATE INDEX idx_ik_case_psss_meddra ON inchikey_case_counts_psss_meddra(inchikey);
            """
        )
        conn.commit()

        print("[INFO] building inchikey_pt_counts_psss_meddra ...", flush=True)
        conn.executescript(
            """
            DROP TABLE IF EXISTS inchikey_pt_counts_psss_meddra;
            CREATE TABLE inchikey_pt_counts_psss_meddra AS
            SELECT
                c.inchikey,
                r.pt_code,
                r.pt_name,
                COUNT(*) AS n_cases
            FROM s4.case_inchikey_psss c
            INNER JOIN reaction_base_meddra r ON r.primaryid = c.primaryid
            GROUP BY c.inchikey, r.pt_code, r.pt_name;

            CREATE INDEX idx_ik_pt_psss_meddra_key ON inchikey_pt_counts_psss_meddra(inchikey);
            CREATE INDEX idx_ik_pt_psss_meddra_pt ON inchikey_pt_counts_psss_meddra(pt_code);
            """
        )
        conn.commit()

        print("[INFO] building model/index/matrix tables ...", flush=True)
        conn.execute("DROP TABLE IF EXISTS inchikey_pt_psss_meddra_model")
        conn.execute(
            """
            CREATE TABLE inchikey_pt_psss_meddra_model AS
            SELECT inchikey, pt_code, pt_name, n_cases
            FROM inchikey_pt_counts_psss_meddra
            WHERE n_cases >= ?
            """,
            (args.min_pair_cases,),
        )
        conn.commit()

        conn.executescript(
            """
            DROP TABLE IF EXISTS inchikey_index_meddra;
            CREATE TABLE inchikey_index_meddra AS
            SELECT
                ROW_NUMBER() OVER (ORDER BY inchikey) AS inchikey_id,
                inchikey
            FROM (SELECT DISTINCT inchikey FROM inchikey_pt_psss_meddra_model);
            CREATE INDEX idx_ik_idx_meddra_key ON inchikey_index_meddra(inchikey);

            DROP TABLE IF EXISTS pt_index_meddra;
            CREATE TABLE pt_index_meddra AS
            SELECT
                ROW_NUMBER() OVER (ORDER BY pt_code) AS pt_id,
                pt_code,
                MIN(pt_name) AS pt_name
            FROM inchikey_pt_psss_meddra_model
            GROUP BY pt_code;
            CREATE INDEX idx_pt_idx_meddra_code ON pt_index_meddra(pt_code);

            DROP TABLE IF EXISTS inchikey_pt_matrix_psss_meddra;
            CREATE TABLE inchikey_pt_matrix_psss_meddra AS
            SELECT
                i.inchikey_id,
                p.pt_id,
                m.n_cases
            FROM inchikey_pt_psss_meddra_model m
            INNER JOIN inchikey_index_meddra i ON i.inchikey = m.inchikey
            INNER JOIN pt_index_meddra p ON p.pt_code = m.pt_code;
            CREATE INDEX idx_mat_meddra_ik ON inchikey_pt_matrix_psss_meddra(inchikey_id);
            CREATE INDEX idx_mat_meddra_pt ON inchikey_pt_matrix_psss_meddra(pt_id);
            """
        )
        conn.commit()

        print("[INFO] exporting step5 structure artifacts ...", flush=True)
        n_case_counts = export_csv(
            conn,
            "SELECT * FROM inchikey_case_counts_psss_meddra ORDER BY n_cases DESC, inchikey",
            args.output_dir / "inchikey_case_counts_psss_meddra.csv",
        )
        n_pair_all = export_csv(
            conn,
            "SELECT * FROM inchikey_pt_counts_psss_meddra ORDER BY n_cases DESC, inchikey, pt_code",
            args.output_dir / "inchikey_pt_counts_psss_meddra.csv",
        )
        n_pair_top = export_csv(
            conn,
            "SELECT * FROM inchikey_pt_counts_psss_meddra ORDER BY n_cases DESC, inchikey, pt_code LIMIT 20000",
            args.output_dir / "inchikey_pt_counts_psss_meddra_top20000.csv",
        )
        n_pair_model = export_csv(
            conn,
            "SELECT * FROM inchikey_pt_psss_meddra_model ORDER BY n_cases DESC, inchikey, pt_code",
            args.output_dir / f"inchikey_pt_counts_psss_meddra_min{args.min_pair_cases}.csv",
        )
        n_ik_idx = export_csv(
            conn,
            "SELECT * FROM inchikey_index_meddra ORDER BY inchikey_id",
            args.output_dir / "inchikey_index_meddra.csv",
        )
        n_pt_idx = export_csv(
            conn,
            "SELECT * FROM pt_index_meddra ORDER BY pt_id",
            args.output_dir / "pt_index_meddra.csv",
        )
        n_edges = export_csv(
            conn,
            "SELECT * FROM inchikey_pt_matrix_psss_meddra ORDER BY inchikey_id, pt_id",
            args.output_dir / f"inchikey_pt_matrix_psss_meddra_min{args.min_pair_cases}.csv",
        )

        finished = dt.datetime.now().isoformat(timespec="seconds")
        report = {
            "started_at": started,
            "finished_at": finished,
            "step4_db": str(args.step4_db),
            "step5_db": str(args.step5_db),
            "min_pair_cases": args.min_pair_cases,
            "case_inchikey_psss_rows": scalar(conn, "SELECT COUNT(*) FROM s4.case_inchikey_psss"),
            "reaction_base_meddra_rows": scalar(conn, "SELECT COUNT(*) FROM reaction_base_meddra"),
            "inchikey_case_counts_psss_meddra_rows": n_case_counts,
            "inchikey_pt_counts_psss_meddra_rows": n_pair_all,
            "inchikey_pt_counts_psss_meddra_top20000_rows": n_pair_top,
            "inchikey_pt_psss_meddra_model_rows": n_pair_model,
            "inchikey_index_meddra_rows": n_ik_idx,
            "pt_index_meddra_rows": n_pt_idx,
            "matrix_edges_rows": n_edges,
            "max_pair_cases_psss_meddra": scalar(
                conn, "SELECT COALESCE(MAX(n_cases),0) FROM inchikey_pt_counts_psss_meddra"
            ),
        }
        (args.output_dir / "step5_structure_report.json").write_text(
            json.dumps(report, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

        conn.execute("PRAGMA wal_checkpoint(TRUNCATE);")
        print("[DONE] MedDRA-standardized InChIKey-PT build completed", flush=True)
        print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    raise SystemExit(main())
