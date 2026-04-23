#!/usr/bin/env python3
"""Step4: build InChI-reaction datasets for modeling (without MedDRA dictionary).

Inputs:
- step2 sqlite: case_drug, case_reaction
- step3 sqlite: drug_term_final (InChI mapping)

Outputs:
- step4 sqlite with aggregated tables
- model-ready sparse matrix CSVs from PS/SS subset
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import sqlite3
from pathlib import Path


def init_db(conn: sqlite3.Connection, rebuild: bool) -> None:
    cur = conn.cursor()
    conn.execute("PRAGMA journal_mode = WAL;")
    conn.execute("PRAGMA synchronous = NORMAL;")
    # Large GROUP BY stages can OOM with MEMORY temp store.
    conn.execute("PRAGMA temp_store = FILE;")

    if rebuild:
        cur.executescript(
            """
            DROP TABLE IF EXISTS term_inchi_map;
            DROP TABLE IF EXISTS inchi_meta;
            DROP TABLE IF EXISTS reaction_base;
            DROP TABLE IF EXISTS case_inchi_any;
            DROP TABLE IF EXISTS case_inchi_psss;
            DROP TABLE IF EXISTS inchi_case_counts_any;
            DROP TABLE IF EXISTS inchi_case_counts_psss;
            DROP TABLE IF EXISTS inchi_pt_counts_any;
            DROP TABLE IF EXISTS inchi_pt_counts_psss;
            DROP TABLE IF EXISTS inchi_pt_psss_model;
            DROP TABLE IF EXISTS inchi_index;
            DROP TABLE IF EXISTS pt_index;
            DROP TABLE IF EXISTS inchi_pt_matrix_psss;
            """
        )
    conn.commit()


def table_exists(conn: sqlite3.Connection, table_name: str) -> bool:
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
        (table_name,),
    ).fetchone()
    return row is not None


def run_sql_stage(
    conn: sqlite3.Connection,
    stage_name: str,
    sql: str,
    produced_tables: list[str],
    resume: bool,
) -> None:
    if resume and all(table_exists(conn, t) for t in produced_tables):
        print(
            f"[SKIP] {stage_name} (existing tables: {', '.join(produced_tables)})",
            flush=True,
        )
        return

    print(f"[INFO] {stage_name} ...", flush=True)
    conn.executescript(sql)
    conn.commit()


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


def scalar(conn: sqlite3.Connection, q: str) -> int:
    return int(conn.execute(q).fetchone()[0])


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build InChI-reaction modeling dataset")
    p.add_argument(
        "--step2-db",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "outputs" / "step2" / "faers_step2.sqlite",
    )
    p.add_argument(
        "--step3-db",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "outputs" / "step3" / "faers_step3.sqlite",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "outputs" / "step4",
    )
    p.add_argument(
        "--min-pair-cases",
        type=int,
        default=20,
        help="Minimum case count for PS/SS InChI-PT pairs in model table",
    )
    p.add_argument("--rebuild", action="store_true")
    p.add_argument(
        "--resume",
        action="store_true",
        help="Resume from existing Step4 tables when possible (skip completed stages)",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if not args.step2_db.exists():
        print(f"[ERROR] step2 db not found: {args.step2_db}")
        return 1
    if not args.step3_db.exists():
        print(f"[ERROR] step3 db not found: {args.step3_db}")
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    step4_db = args.output_dir / "faers_step4.sqlite"

    started = dt.datetime.now().isoformat(timespec="seconds")

    conn = sqlite3.connect(step4_db)
    try:
        init_db(conn, rebuild=args.rebuild)

        conn.execute("ATTACH DATABASE ? AS s2", (str(args.step2_db),))
        conn.execute("ATTACH DATABASE ? AS s3", (str(args.step3_db),))

        if args.rebuild and args.resume:
            print("[WARN] --rebuild is set, --resume is ignored.", flush=True)
        resume_mode = args.resume and (not args.rebuild)

        run_sql_stage(
            conn,
            "building term_inchi_map",
            """
            DROP TABLE IF EXISTS term_inchi_map;
            CREATE TABLE term_inchi_map AS
            SELECT DISTINCT
                term_std,
                inchi,
                dbid,
                drugbank_name,
                rxcui
            FROM s3.drug_term_final
            WHERE TRIM(COALESCE(inchi,'')) <> '';

            CREATE INDEX idx_term_inchi_term ON term_inchi_map(term_std);
            CREATE INDEX idx_term_inchi_inchi ON term_inchi_map(inchi);
            """,
            produced_tables=["term_inchi_map"],
            resume=resume_mode,
        )

        run_sql_stage(
            conn,
            "building inchi metadata",
            """
            DROP TABLE IF EXISTS inchi_meta;
            CREATE TABLE inchi_meta AS
            SELECT
                inchi,
                MIN(dbid) AS dbid_any,
                MIN(drugbank_name) AS drugbank_name_any,
                COUNT(DISTINCT dbid) AS dbid_count,
                COUNT(DISTINCT term_std) AS mapped_term_count
            FROM term_inchi_map
            GROUP BY inchi;

            CREATE INDEX idx_inchi_meta_inchi ON inchi_meta(inchi);
            """,
            produced_tables=["inchi_meta"],
            resume=resume_mode,
        )

        run_sql_stage(
            conn,
            "preparing reaction_base",
            """
            DROP TABLE IF EXISTS reaction_base;
            CREATE TABLE reaction_base AS
            SELECT DISTINCT
                primaryid,
                caseid,
                pt_norm
            FROM s2.case_reaction
            WHERE TRIM(COALESCE(pt_norm,'')) <> '';

            CREATE INDEX idx_reaction_base_primaryid ON reaction_base(primaryid);
            CREATE INDEX idx_reaction_base_pt ON reaction_base(pt_norm);
            """,
            produced_tables=["reaction_base"],
            resume=resume_mode,
        )

        run_sql_stage(
            conn,
            "building case_inchi_any / case_inchi_psss",
            """
            DROP TABLE IF EXISTS case_inchi_any;
            CREATE TABLE case_inchi_any AS
            SELECT DISTINCT
                d.primaryid,
                d.caseid,
                m.inchi
            FROM s2.case_drug d
            INNER JOIN term_inchi_map m
                ON d.drug_key_std = m.term_std
            WHERE TRIM(COALESCE(d.primaryid,'')) <> '';

            CREATE INDEX idx_case_inchi_any_primaryid ON case_inchi_any(primaryid);
            CREATE INDEX idx_case_inchi_any_inchi ON case_inchi_any(inchi);

            DROP TABLE IF EXISTS case_inchi_psss;
            CREATE TABLE case_inchi_psss AS
            SELECT DISTINCT
                d.primaryid,
                d.caseid,
                m.inchi
            FROM s2.case_drug d
            INNER JOIN term_inchi_map m
                ON d.drug_key_std = m.term_std
            WHERE TRIM(COALESCE(d.primaryid,'')) <> ''
              AND UPPER(TRIM(COALESCE(d.role_cod,''))) IN ('PS', 'SS');

            CREATE INDEX idx_case_inchi_psss_primaryid ON case_inchi_psss(primaryid);
            CREATE INDEX idx_case_inchi_psss_inchi ON case_inchi_psss(inchi);
            """,
            produced_tables=["case_inchi_any", "case_inchi_psss"],
            resume=resume_mode,
        )

        run_sql_stage(
            conn,
            "aggregating inchi-level case counts",
            """
            DROP TABLE IF EXISTS inchi_case_counts_any;
            CREATE TABLE inchi_case_counts_any AS
            SELECT
                inchi,
                COUNT(*) AS n_cases
            FROM case_inchi_any
            GROUP BY inchi;

            CREATE INDEX idx_inchi_case_any_inchi ON inchi_case_counts_any(inchi);

            DROP TABLE IF EXISTS inchi_case_counts_psss;
            CREATE TABLE inchi_case_counts_psss AS
            SELECT
                inchi,
                COUNT(*) AS n_cases
            FROM case_inchi_psss
            GROUP BY inchi;

            CREATE INDEX idx_inchi_case_psss_inchi ON inchi_case_counts_psss(inchi);
            """,
            produced_tables=["inchi_case_counts_any", "inchi_case_counts_psss"],
            resume=resume_mode,
        )

        run_sql_stage(
            conn,
            "aggregating InChI-PT counts (ANY)",
            """
            DROP TABLE IF EXISTS inchi_pt_counts_any;
            CREATE TABLE inchi_pt_counts_any AS
            SELECT
                c.inchi,
                r.pt_norm,
                COUNT(*) AS n_cases
            FROM case_inchi_any c
            INNER JOIN reaction_base r
                ON r.primaryid = c.primaryid
            GROUP BY c.inchi, r.pt_norm;

            CREATE INDEX idx_inchi_pt_any_inchi ON inchi_pt_counts_any(inchi);
            CREATE INDEX idx_inchi_pt_any_pt ON inchi_pt_counts_any(pt_norm);
            """,
            produced_tables=["inchi_pt_counts_any"],
            resume=resume_mode,
        )

        run_sql_stage(
            conn,
            "aggregating InChI-PT counts (PS/SS)",
            """
            DROP TABLE IF EXISTS inchi_pt_counts_psss;
            CREATE TABLE inchi_pt_counts_psss AS
            SELECT
                c.inchi,
                r.pt_norm,
                COUNT(*) AS n_cases
            FROM case_inchi_psss c
            INNER JOIN reaction_base r
                ON r.primaryid = c.primaryid
            GROUP BY c.inchi, r.pt_norm;

            CREATE INDEX idx_inchi_pt_psss_inchi ON inchi_pt_counts_psss(inchi);
            CREATE INDEX idx_inchi_pt_psss_pt ON inchi_pt_counts_psss(pt_norm);
            """,
            produced_tables=["inchi_pt_counts_psss"],
            resume=resume_mode,
        )

        rebuild_model_tables = not (
            resume_mode
            and table_exists(conn, "inchi_pt_psss_model")
            and table_exists(conn, "inchi_index")
            and table_exists(conn, "pt_index")
            and table_exists(conn, "inchi_pt_matrix_psss")
        )
        if rebuild_model_tables:
            print("[INFO] creating model table and sparse matrix ...", flush=True)
            conn.execute("DROP TABLE IF EXISTS inchi_pt_psss_model")
            conn.execute(
                """
                CREATE TABLE inchi_pt_psss_model AS
                SELECT inchi, pt_norm, n_cases
                FROM inchi_pt_counts_psss
                WHERE n_cases >= ?
                """,
                (args.min_pair_cases,),
            )
            conn.commit()
        else:
            print("[SKIP] creating model table and sparse matrix (already exists)", flush=True)

        if rebuild_model_tables:
            conn.executescript(
                """
                DROP TABLE IF EXISTS inchi_index;
                CREATE TABLE inchi_index AS
                SELECT
                    ROW_NUMBER() OVER (ORDER BY inchi) AS inchi_id,
                    inchi
                FROM (SELECT DISTINCT inchi FROM inchi_pt_psss_model);
                CREATE INDEX idx_inchi_index_inchi ON inchi_index(inchi);

                DROP TABLE IF EXISTS pt_index;
                CREATE TABLE pt_index AS
                SELECT
                    ROW_NUMBER() OVER (ORDER BY pt_norm) AS pt_id,
                    pt_norm
                FROM (SELECT DISTINCT pt_norm FROM inchi_pt_psss_model);
                CREATE INDEX idx_pt_index_pt ON pt_index(pt_norm);

                DROP TABLE IF EXISTS inchi_pt_matrix_psss;
                CREATE TABLE inchi_pt_matrix_psss AS
                SELECT
                    i.inchi_id,
                    p.pt_id,
                    m.n_cases
                FROM inchi_pt_psss_model m
                INNER JOIN inchi_index i ON i.inchi = m.inchi
                INNER JOIN pt_index p ON p.pt_norm = m.pt_norm;

                CREATE INDEX idx_matrix_inchi ON inchi_pt_matrix_psss(inchi_id);
                CREATE INDEX idx_matrix_pt ON inchi_pt_matrix_psss(pt_id);
                """
            )
            conn.commit()

        # Export
        n_meta = export_csv(
            conn,
            "SELECT * FROM inchi_meta ORDER BY mapped_term_count DESC, inchi",
            args.output_dir / "inchi_meta.csv",
        )
        n_inchi_any = export_csv(
            conn,
            "SELECT * FROM inchi_case_counts_any ORDER BY n_cases DESC, inchi",
            args.output_dir / "inchi_case_counts_any.csv",
        )
        n_inchi_psss = export_csv(
            conn,
            "SELECT * FROM inchi_case_counts_psss ORDER BY n_cases DESC, inchi",
            args.output_dir / "inchi_case_counts_psss.csv",
        )
        n_pair_psss = export_csv(
            conn,
            "SELECT * FROM inchi_pt_psss_model ORDER BY n_cases DESC, inchi, pt_norm",
            args.output_dir / f"inchi_pt_counts_psss_min{args.min_pair_cases}.csv",
        )
        n_inchi_idx = export_csv(
            conn,
            "SELECT * FROM inchi_index ORDER BY inchi_id",
            args.output_dir / "inchi_index.csv",
        )
        n_pt_idx = export_csv(
            conn,
            "SELECT * FROM pt_index ORDER BY pt_id",
            args.output_dir / "pt_index.csv",
        )
        n_edges = export_csv(
            conn,
            "SELECT * FROM inchi_pt_matrix_psss ORDER BY inchi_id, pt_id",
            args.output_dir / f"inchi_pt_matrix_psss_min{args.min_pair_cases}.csv",
        )

        # Additional top tables for quick review
        export_csv(
            conn,
            "SELECT * FROM inchi_pt_counts_any ORDER BY n_cases DESC, inchi, pt_norm LIMIT 20000",
            args.output_dir / "inchi_pt_counts_any_top20000.csv",
        )
        export_csv(
            conn,
            "SELECT * FROM inchi_pt_counts_psss ORDER BY n_cases DESC, inchi, pt_norm LIMIT 20000",
            args.output_dir / "inchi_pt_counts_psss_top20000.csv",
        )

        finished = dt.datetime.now().isoformat(timespec="seconds")
        summary = {
            "started_at": started,
            "finished_at": finished,
            "min_pair_cases": args.min_pair_cases,
            "term_inchi_map_rows": scalar(conn, "SELECT COUNT(*) FROM term_inchi_map"),
            "inchi_meta_rows": n_meta,
            "reaction_base_rows": scalar(conn, "SELECT COUNT(*) FROM reaction_base"),
            "case_inchi_any_rows": scalar(conn, "SELECT COUNT(*) FROM case_inchi_any"),
            "case_inchi_psss_rows": scalar(conn, "SELECT COUNT(*) FROM case_inchi_psss"),
            "inchi_case_any_rows": n_inchi_any,
            "inchi_case_psss_rows": n_inchi_psss,
            "inchi_pt_any_rows": scalar(conn, "SELECT COUNT(*) FROM inchi_pt_counts_any"),
            "inchi_pt_psss_rows": scalar(conn, "SELECT COUNT(*) FROM inchi_pt_counts_psss"),
            "inchi_pt_psss_model_rows": n_pair_psss,
            "inchi_index_rows": n_inchi_idx,
            "pt_index_rows": n_pt_idx,
            "matrix_edges_rows": n_edges,
            "max_pair_cases_psss": scalar(conn, "SELECT COALESCE(MAX(n_cases),0) FROM inchi_pt_counts_psss"),
            "max_pair_cases_any": scalar(conn, "SELECT COALESCE(MAX(n_cases),0) FROM inchi_pt_counts_any"),
        }

        (args.output_dir / "step4_report.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

        # Merge WAL back into main DB so downstream tools see a stable file size.
        conn.execute("PRAGMA wal_checkpoint(TRUNCATE);")

        print("[DONE] Step4 completed", flush=True)
        print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
        return 0
    except Exception as exc:
        failed = dt.datetime.now().isoformat(timespec="seconds")
        error_summary = {
            "started_at": started,
            "failed_at": failed,
            "error": str(exc),
        }
        (args.output_dir / "step4_error.json").write_text(
            json.dumps(error_summary, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        print("[ERROR] Step4 failed", flush=True)
        print(json.dumps(error_summary, ensure_ascii=False, indent=2), flush=True)
        raise
    finally:
        conn.close()


if __name__ == "__main__":
    raise SystemExit(main())

