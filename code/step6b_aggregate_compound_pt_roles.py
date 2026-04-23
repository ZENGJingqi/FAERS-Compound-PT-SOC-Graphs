from __future__ import annotations

import json
import sqlite3
import time
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "outputs" / "step6_core_graphs" / "core_graph_full"
OUT_DB = OUT_DIR / "step6_compound_pt_soc_core_full.sqlite"
OUT_DB_REL = "outputs/step6_core_graphs/core_graph_full/step6_compound_pt_soc_core_full.sqlite"
STEP5_DB = ROOT / "outputs" / "step5" / "faers_step5.sqlite"
STEP4_DB = ROOT / "outputs" / "step4" / "faers_step4.sqlite"
STEP5_DB_REL = "outputs/step5/faers_step5.sqlite"
MEDDRA_MDHIER = ROOT.parents[1] / "reference_data" / "MedDRA" / "MedDRA_29_0_English" / "MedAscii" / "mdhier.asc"

BATCH_SIZE = 100_000


def attach_path(path: Path) -> str:
    return str(path).replace("\\", "/")


def scalar(conn: sqlite3.Connection, q: str):
    return conn.execute(q).fetchone()[0]


def load_pt_soc_from_meddra() -> pd.DataFrame:
    rows_primary = []
    rows_fallback = []
    with MEDDRA_MDHIER.open("r", encoding="latin-1") as f:
        for line in f:
            parts = line.rstrip("\n\r").split("$")
            if len(parts) < 12:
                continue
            pt_code = parts[0].strip()
            pt_name = parts[4].strip()
            soc_code = parts[3].strip()
            soc_name = parts[7].strip()
            primary_soc_code = parts[10].strip() or soc_code
            primary_flag = parts[11].strip().upper()
            row = {
                "pt_code": pt_code,
                "pt_name": pt_name,
                "soc_code_primary": primary_soc_code,
                "soc_name_primary": soc_name,
            }
            rows_fallback.append(row)
            if primary_flag == "Y":
                rows_primary.append(row)

    primary_df = pd.DataFrame(rows_primary).drop_duplicates("pt_code")
    fallback_df = pd.DataFrame(rows_fallback).drop_duplicates("pt_code")
    missing = fallback_df.loc[~fallback_df["pt_code"].isin(primary_df["pt_code"])]
    pt_soc_df = pd.concat([primary_df, missing], ignore_index=True).drop_duplicates("pt_code")
    return pt_soc_df


def ensure_node_tables(conn: sqlite3.Connection) -> None:
    compound_df = pd.read_sql_query(
        """
        SELECT *
        FROM s4.inchikey_meta
        WHERE inchikey IN (
            SELECT inchikey FROM report_compound_ps
            UNION
            SELECT inchikey FROM report_compound_ss
        )
        """,
        conn,
    )
    pt_conn = sqlite3.connect(STEP5_DB_REL)
    pt_df = pd.read_sql_query(
        """
        SELECT pt_code, pt_name
        FROM pt_index_meddra
        """,
        pt_conn,
    )
    pt_conn.close()

    pt_soc_df = load_pt_soc_from_meddra()
    pt_codes_in_data = set(pt_df["pt_code"].astype(str))
    pt_soc_df["pt_code"] = pt_soc_df["pt_code"].astype(str)
    pt_soc_df = pt_soc_df[pt_soc_df["pt_code"].isin(pt_codes_in_data)].reset_index(drop=True)
    soc_df = (
        pt_soc_df[["soc_name_primary", "soc_code_primary"]]
        .dropna(subset=["soc_name_primary"])
        .drop_duplicates()
        .rename(columns={"soc_name_primary": "soc_name", "soc_code_primary": "soc_code"})
        .reset_index(drop=True)
    )

    for table in ["node_compound", "node_pt", "node_soc", "edge_pt_soc_primary"]:
        conn.execute(f"DROP TABLE IF EXISTS {table}")

    compound_df.to_sql("node_compound", conn, index=False, if_exists="replace")
    pt_df.to_sql("node_pt", conn, index=False, if_exists="replace")
    soc_df.to_sql("node_soc", conn, index=False, if_exists="replace")
    pt_soc_df.to_sql("edge_pt_soc_primary", conn, index=False, if_exists="replace")

    conn.execute("CREATE UNIQUE INDEX idx_node_compound_key ON node_compound(inchikey)")
    conn.execute("CREATE UNIQUE INDEX idx_node_pt_code ON node_pt(pt_code)")
    conn.execute("CREATE UNIQUE INDEX idx_node_soc_name ON node_soc(soc_name)")
    conn.execute("CREATE INDEX idx_edge_pt_soc_pt ON edge_pt_soc_primary(pt_code)")
    conn.execute("CREATE INDEX idx_edge_pt_soc_soc ON edge_pt_soc_primary(soc_name_primary)")
    conn.commit()


def aggregate_role(conn: sqlite3.Connection, role: str) -> dict:
    source_table = f"report_compound_{role.lower()}"
    target_table = f"edge_compound_pt_{role.lower()}"
    conn.execute(f"DROP TABLE IF EXISTS {target_table}")
    conn.execute(
        f"""
        CREATE TABLE {target_table} (
            inchikey TEXT NOT NULL,
            pt_code TEXT NOT NULL,
            pt_name TEXT,
            n_reports INTEGER NOT NULL,
            PRIMARY KEY(inchikey, pt_code)
        ) WITHOUT ROWID
        """
    )
    conn.commit()

    read_conn = sqlite3.connect(OUT_DB_REL)
    read_conn.execute("PRAGMA journal_mode=WAL;")
    total_source = scalar(read_conn, f"SELECT COUNT(*) FROM {source_table}")
    cur = read_conn.cursor()
    cur.execute(f"SELECT primaryid, inchikey FROM {source_table} ORDER BY primaryid")
    processed = 0
    batch_no = 0
    t0 = time.time()

    while True:
        rows = cur.fetchmany(BATCH_SIZE)
        if not rows:
            break
        batch_no += 1
        processed += len(rows)

        conn.execute("DROP TABLE IF EXISTS temp.batch_report_compound")
        conn.execute("CREATE TEMP TABLE batch_report_compound(primaryid TEXT, inchikey TEXT)")
        conn.executemany("INSERT INTO batch_report_compound(primaryid, inchikey) VALUES (?, ?)", rows)
        conn.execute("CREATE INDEX idx_batch_report_compound_pid ON batch_report_compound(primaryid)")

        conn.execute(
            f"""
            INSERT INTO {target_table}(inchikey, pt_code, pt_name, n_reports)
            SELECT
                b.inchikey,
                r.pt_code,
                MIN(r.pt_name) AS pt_name,
                COUNT(*) AS n_reports
            FROM batch_report_compound b
            INNER JOIN s5.reaction_base_meddra r
                ON b.primaryid = r.primaryid
            GROUP BY b.inchikey, r.pt_code
            ON CONFLICT(inchikey, pt_code) DO UPDATE SET
                n_reports = {target_table}.n_reports + excluded.n_reports,
                pt_name = COALESCE({target_table}.pt_name, excluded.pt_name)
            """
        )
        conn.commit()
        if batch_no % 5 == 0:
            current_edges = scalar(conn, f"SELECT COUNT(*) FROM {target_table}")
            print(
                f"[INFO] {role}: processed {processed}/{total_source} report-compound rows; current compound-pt edges={current_edges}; elapsed={round(time.time()-t0,1)}s",
                flush=True,
            )

    conn.execute(f"CREATE INDEX idx_{target_table}_inchikey ON {target_table}(inchikey)")
    conn.execute(f"CREATE INDEX idx_{target_table}_ptcode ON {target_table}(pt_code)")
    conn.commit()
    read_conn.close()

    return {
        "role": role,
        "source_rows": int(total_source),
        "edge_rows": int(scalar(conn, f"SELECT COUNT(*) FROM {target_table}")),
        "total_weight": int(scalar(conn, f"SELECT SUM(n_reports) FROM {target_table}")),
        "unique_compounds": int(scalar(conn, f"SELECT COUNT(DISTINCT inchikey) FROM {target_table}")),
        "unique_pts": int(scalar(conn, f"SELECT COUNT(DISTINCT pt_code) FROM {target_table}")),
        "elapsed_sec": round(time.time() - t0, 1),
    }


def main() -> None:
    conn = sqlite3.connect(OUT_DB)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=OFF;")
    conn.execute("PRAGMA temp_store=MEMORY;")
    conn.execute("PRAGMA cache_size=-200000;")
    conn.execute(f"ATTACH DATABASE '{attach_path(STEP5_DB)}' AS s5;")
    conn.execute(f"ATTACH DATABASE '{attach_path(STEP4_DB)}' AS s4;")

    print("[INFO] ensuring node/PT/SOC tables ...", flush=True)
    ensure_node_tables(conn)

    ps = aggregate_role(conn, "PS")
    ss = aggregate_role(conn, "SS")

    node_summary = pd.DataFrame(
        [
            ("compound", scalar(conn, "SELECT COUNT(*) FROM node_compound")),
            ("pt", scalar(conn, "SELECT COUNT(*) FROM node_pt")),
            ("soc", scalar(conn, "SELECT COUNT(*) FROM node_soc")),
        ],
        columns=["node_type", "n_nodes"],
    )
    edge_summary = pd.DataFrame(
        [
            ("compound_has_pt_ps", "compound", "pt", scalar(conn, "SELECT COUNT(*) FROM edge_compound_pt_ps")),
            ("compound_has_pt_ss", "compound", "pt", scalar(conn, "SELECT COUNT(*) FROM edge_compound_pt_ss")),
            ("pt_belongs_to_primary_soc", "pt", "soc", scalar(conn, "SELECT COUNT(*) FROM edge_pt_soc_primary")),
        ],
        columns=["edge_type", "src_type", "dst_type", "n_edges"],
    )
    role_summary = pd.DataFrame(
        [
            ("ps_total_weight", ps["total_weight"]),
            ("ss_total_weight", ss["total_weight"]),
            ("ps_unique_reports", scalar(conn, "SELECT COUNT(DISTINCT primaryid) FROM report_compound_ps")),
            ("ss_unique_reports", scalar(conn, "SELECT COUNT(DISTINCT primaryid) FROM report_compound_ss")),
        ],
        columns=["metric", "value"],
    )
    node_summary.to_csv(OUT_DIR / "node_summary.csv", index=False)
    edge_summary.to_csv(OUT_DIR / "edge_summary.csv", index=False)
    role_summary.to_csv(OUT_DIR / "role_summary.csv", index=False)

    summary = {"PS": ps, "SS": ss}
    (OUT_DIR / "compound_pt_role_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    conn.close()
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
