from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
CORE_DB = ROOT / "outputs" / "step14_heterograph" / "compound_pt_role_core_v3" / "step14_compound_pt_role_core_v3.sqlite"
TEMP_DB = ROOT / "outputs" / "step15_temporal_link_dataset" / "step15_temporal_link_dataset.sqlite"
MEDDRA_MDHIER = ROOT.parent / "MedDRA" / "MedDRA_29_0_English" / "MedAscii" / "mdhier.asc"


def load_pt_soc() -> pd.DataFrame:
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
    return pd.concat([primary_df, missing], ignore_index=True).drop_duplicates("pt_code")


def repair_db(db_path: Path, edge_tables: list[str]) -> dict:
    conn = sqlite3.connect(db_path)
    union_sql = " UNION ".join([f"SELECT pt_code, pt_name FROM {t}" for t in edge_tables])
    pt_df = pd.read_sql_query(
        f"SELECT pt_code, MIN(pt_name) AS pt_name FROM ({union_sql}) GROUP BY pt_code",
        conn,
    )
    pt_df["pt_code"] = pt_df["pt_code"].astype(str)
    pt_soc_df = load_pt_soc()
    pt_soc_df["pt_code"] = pt_soc_df["pt_code"].astype(str)
    pt_soc_df = pt_soc_df[pt_soc_df["pt_code"].isin(set(pt_df["pt_code"]))].reset_index(drop=True)
    soc_df = (
        pt_soc_df[["soc_name_primary", "soc_code_primary"]]
        .drop_duplicates()
        .rename(columns={"soc_name_primary": "soc_name", "soc_code_primary": "soc_code"})
        .reset_index(drop=True)
    )

    for t in ["node_pt", "edge_pt_soc_primary", "node_soc"]:
        conn.execute(f"DROP TABLE IF EXISTS {t}")
    pt_df.to_sql("node_pt", conn, index=False, if_exists="replace")
    pt_soc_df[["pt_code", "soc_code_primary", "soc_name_primary"]].to_sql(
        "edge_pt_soc_primary", conn, index=False, if_exists="replace"
    )
    soc_df.to_sql("node_soc", conn, index=False, if_exists="replace")
    conn.execute("CREATE UNIQUE INDEX IF NOT EXISTS idx_node_pt_code ON node_pt(pt_code)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_edge_pt_soc_pt ON edge_pt_soc_primary(pt_code)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_edge_pt_soc_soc ON edge_pt_soc_primary(soc_name_primary)")
    conn.execute("CREATE UNIQUE INDEX IF NOT EXISTS idx_node_soc_name ON node_soc(soc_name)")
    conn.commit()
    summary = {
        "db": str(db_path),
        "node_pt": int(conn.execute("select count(*) from node_pt").fetchone()[0]),
        "edge_pt_soc_primary": int(conn.execute("select count(*) from edge_pt_soc_primary").fetchone()[0]),
        "node_soc": int(conn.execute("select count(*) from node_soc").fetchone()[0]),
    }
    conn.close()
    return summary


def write_summaries():
    core_conn = sqlite3.connect(CORE_DB)
    core_node_summary = pd.DataFrame(
        [
            ("compound", core_conn.execute("select count(*) from node_compound").fetchone()[0]),
            ("pt", core_conn.execute("select count(*) from node_pt").fetchone()[0]),
            ("soc", core_conn.execute("select count(*) from node_soc").fetchone()[0]),
        ],
        columns=["node_type", "n_nodes"],
    )
    core_edge_summary = pd.DataFrame(
        [
            ("compound_has_pt_ps", "compound", "pt", core_conn.execute("select count(*) from edge_compound_pt_ps").fetchone()[0]),
            ("compound_has_pt_ss", "compound", "pt", core_conn.execute("select count(*) from edge_compound_pt_ss").fetchone()[0]),
            ("pt_belongs_to_primary_soc", "pt", "soc", core_conn.execute("select count(*) from edge_pt_soc_primary").fetchone()[0]),
        ],
        columns=["edge_type", "src_type", "dst_type", "n_edges"],
    )
    core_node_summary.to_csv(CORE_DB.parent / "node_summary.csv", index=False)
    core_edge_summary.to_csv(CORE_DB.parent / "edge_summary.csv", index=False)
    core_conn.close()


def main():
    core = repair_db(CORE_DB, ["edge_compound_pt_ps", "edge_compound_pt_ss"])
    temp = repair_db(
        TEMP_DB,
        ["edge_compound_pt_ps_train", "edge_compound_pt_ps_val", "edge_compound_pt_ps_test", "edge_compound_pt_ss_train", "edge_compound_pt_ss_val", "edge_compound_pt_ss_test"],
    )
    write_summaries()
    print({"core": core, "temporal": temp})


if __name__ == "__main__":
    main()
