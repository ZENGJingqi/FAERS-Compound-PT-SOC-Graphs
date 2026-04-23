from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
MEDDRA_ZH = ROOT.parents[1] / "reference_data" / "MedDRA" / "MedDRA_29_0_Chinese" / "ascii-290"
MEDDRA_EN = ROOT.parents[1] / "reference_data" / "MedDRA" / "MedDRA_29_0_English" / "MedAscii"
STEP6_DIR = ROOT / "outputs" / "step6_core_graphs"
BASE_DB = STEP6_DIR / "core_graph_full" / "step6_compound_pt_soc_core_full.sqlite"
GRAPH_DIRS = {
    "ge10": STEP6_DIR / "core_graph_ge10" / "step6_compound_pt_soc_core_ge10.sqlite",
    "ge20": STEP6_DIR / "core_graph_ge20" / "step6_compound_pt_soc_core_ge20.sqlite",
    "ge30": STEP6_DIR / "core_graph_ge30" / "step6_compound_pt_soc_core_ge30.sqlite",
}
REFERENCE_DIR = ROOT / "outputs" / "reference_tables"


def read_asc_map(path: Path, key_idx: int = 0, val_idx: int = 1) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.rstrip("\r\n").split("$")
            if len(parts) <= max(key_idx, val_idx):
                continue
            key = parts[key_idx].strip()
            val = parts[val_idx].strip()
            if key:
                out[key] = val
    return out


def ensure_column(conn: sqlite3.Connection, table: str, col: str, col_type: str) -> None:
    cols = {row[1] for row in conn.execute(f"PRAGMA table_info({table})")}
    if col not in cols:
        conn.execute(f"ALTER TABLE {table} ADD COLUMN {col} {col_type}")


def fill_graph_bilingual(db_path: Path, pt_en: dict[str, str], pt_zh: dict[str, str], soc_en: dict[str, str], soc_zh: dict[str, str]) -> None:
    conn = sqlite3.connect(str(db_path))
    try:
        table_columns = {
            "node_pt": [("pt_name_en", "TEXT"), ("pt_name_zh", "TEXT")],
            "edge_compound_pt_ps": [("pt_name_en", "TEXT"), ("pt_name_zh", "TEXT")],
            "edge_compound_pt_ss": [("pt_name_en", "TEXT"), ("pt_name_zh", "TEXT")],
            "node_soc": [("soc_name_en", "TEXT"), ("soc_name_zh", "TEXT")],
            "edge_pt_soc_primary": [("soc_name_en", "TEXT"), ("soc_name_zh", "TEXT")],
        }
        for table, cols in table_columns.items():
            for col, typ in cols:
                ensure_column(conn, table, col, typ)

        node_pt_rows = [(pt_en.get(code, name), pt_zh.get(code, ""), code) for code, name in conn.execute("SELECT pt_code, pt_name FROM node_pt")]
        conn.executemany("UPDATE node_pt SET pt_name_en=?, pt_name_zh=? WHERE pt_code=?", node_pt_rows)

        for table in ["edge_compound_pt_ps", "edge_compound_pt_ss"]:
            edge_rows = [(pt_en.get(code, name), pt_zh.get(code, ""), code) for code, name in conn.execute(f"SELECT DISTINCT pt_code, pt_name FROM {table}")]
            conn.executemany(f"UPDATE {table} SET pt_name_en=?, pt_name_zh=? WHERE pt_code=?", edge_rows)

        node_soc_rows = [(soc_en.get(code, name), soc_zh.get(code, ""), code) for name, code in conn.execute("SELECT soc_name, soc_code FROM node_soc")]
        conn.executemany("UPDATE node_soc SET soc_name_en=?, soc_name_zh=? WHERE soc_code=?", node_soc_rows)

        edge_soc_rows = [(soc_en.get(code, name), soc_zh.get(code, ""), code) for code, name in conn.execute("SELECT DISTINCT soc_code_primary, soc_name_primary FROM edge_pt_soc_primary")]
        conn.executemany("UPDATE edge_pt_soc_primary SET soc_name_en=?, soc_name_zh=? WHERE soc_code_primary=?", edge_soc_rows)
        conn.commit()
    finally:
        conn.close()


def build_reference_tables(pt_en: dict[str, str], pt_zh: dict[str, str], soc_en: dict[str, str], soc_zh: dict[str, str]) -> None:
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    base = sqlite3.connect(str(BASE_DB))
    try:
        compounds = pd.read_sql_query("SELECT inchikey, smiles_std, inchi_std, term_count, dbid_count, dbid_list, drugbank_name_list FROM node_compound ORDER BY inchikey", base)
        reactions = pd.read_sql_query("SELECT DISTINCT p.pt_code, p.pt_name AS pt_name_base, e.soc_code_primary, e.soc_name_primary FROM node_pt p LEFT JOIN edge_pt_soc_primary e ON p.pt_code = e.pt_code ORDER BY p.pt_code", base)
        ps_assoc = pd.read_sql_query("SELECT inchikey, pt_code, pt_name, n_reports, pt_name_en, pt_name_zh FROM edge_compound_pt_ps ORDER BY n_reports DESC, inchikey, pt_code", base)
        ss_assoc = pd.read_sql_query("SELECT inchikey, pt_code, pt_name, n_reports, pt_name_en, pt_name_zh FROM edge_compound_pt_ss ORDER BY n_reports DESC, inchikey, pt_code", base)
    finally:
        base.close()

    graph_compounds: dict[str, set[str]] = {}
    graph_pts: dict[str, set[str]] = {}
    for name, path in GRAPH_DIRS.items():
        conn = sqlite3.connect(str(path))
        try:
            graph_compounds[name] = {row[0] for row in conn.execute("SELECT inchikey FROM node_compound")}
            graph_pts[name] = {row[0] for row in conn.execute("SELECT pt_code FROM node_pt")}
        finally:
            conn.close()

    for name in ["ge10", "ge20", "ge30"]:
        compounds[f"in_{name}"] = compounds["inchikey"].map(lambda x: 1 if x in graph_compounds[name] else 0)

    reactions["pt_name_en"] = reactions.apply(lambda r: pt_en.get(r["pt_code"], r["pt_name_base"] if pd.notna(r["pt_name_base"]) else ""), axis=1)
    reactions["pt_name_zh"] = reactions["pt_code"].map(lambda x: pt_zh.get(x, ""))
    reactions["soc_name_en"] = reactions["soc_code_primary"].map(lambda x: soc_en.get(x, "") if pd.notna(x) else "")
    reactions["soc_name_zh"] = reactions["soc_code_primary"].map(lambda x: soc_zh.get(x, "") if pd.notna(x) else "")
    for name in ["ge10", "ge20", "ge30"]:
        reactions[f"in_{name}"] = reactions["pt_code"].map(lambda x: 1 if x in graph_pts[name] else 0)

    reactions = reactions[["pt_code", "pt_name_en", "pt_name_zh", "soc_code_primary", "soc_name_en", "soc_name_zh", "in_ge10", "in_ge20", "in_ge30"]]
    soc_df = reactions[["soc_code_primary", "soc_name_en", "soc_name_zh"]].dropna().drop_duplicates().sort_values(["soc_code_primary"])

    with pd.ExcelWriter(REFERENCE_DIR / "all_compounds_basic_info.xlsx", engine="openpyxl") as writer:
        compounds.to_excel(writer, sheet_name="compounds", index=False)
    with pd.ExcelWriter(REFERENCE_DIR / "all_reactions_basic_info_bilingual.xlsx", engine="openpyxl") as writer:
        reactions.to_excel(writer, sheet_name="reactions", index=False)
    with pd.ExcelWriter(REFERENCE_DIR / "all_soc_basic_info_bilingual.xlsx", engine="openpyxl") as writer:
        soc_df.to_excel(writer, sheet_name="soc", index=False)

    ps_assoc.to_csv(REFERENCE_DIR / "all_compound_pt_ps_associations_bilingual.csv.gz", index=False, compression="gzip")
    ss_assoc.to_csv(REFERENCE_DIR / "all_compound_pt_ss_associations_bilingual.csv.gz", index=False, compression="gzip")

    summary = {
        "compound_rows": int(len(compounds)),
        "reaction_rows": int(len(reactions)),
        "soc_rows": int(len(soc_df)),
        "ps_association_rows": int(len(ps_assoc)),
        "ss_association_rows": int(len(ss_assoc)),
    }
    (REFERENCE_DIR / "reference_tables_summary.json").write_text(__import__('json').dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")


def main() -> None:
    pt_en = read_asc_map(MEDDRA_EN / "pt.asc")
    pt_zh = read_asc_map(MEDDRA_ZH / "pt.asc")
    soc_en = read_asc_map(MEDDRA_EN / "soc.asc")
    soc_zh = read_asc_map(MEDDRA_ZH / "soc.asc")
    for db in [BASE_DB, *GRAPH_DIRS.values()]:
        fill_graph_bilingual(db, pt_en, pt_zh, soc_en, soc_zh)
    build_reference_tables(pt_en, pt_zh, soc_en, soc_zh)
    print("done")


if __name__ == "__main__":
    main()

