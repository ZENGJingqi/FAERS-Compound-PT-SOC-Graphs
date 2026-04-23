from __future__ import annotations

import shutil
import sqlite3
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
MEDDRA_ZH = ROOT.parent / "MedDRA" / "MedDRA_29_0_Chinese" / "ascii-290"
MEDDRA_EN = ROOT.parent / "MedDRA" / "MedDRA_29_0_English" / "MedAscii"
STEP14_DIR = ROOT / "outputs" / "step14_heterograph"
BASE_DB = STEP14_DIR / "compound_pt_role_core_v3" / "step14_compound_pt_role_core_v3.sqlite"
GRAPH_DIRS = {
    "ge10": STEP14_DIR / "compound_pt_role_core_ge10" / "step14_compound_pt_role_core_ge10.sqlite",
    "ge20": STEP14_DIR / "compound_pt_role_core_ge20" / "step14_compound_pt_role_core_ge20.sqlite",
    "ge30": STEP14_DIR / "compound_pt_role_core_ge30" / "step14_compound_pt_role_core_ge30.sqlite",
}
PACKAGE_DIR = ROOT / "expert_package_cleaning_and_core_graphs_2026-04-22"


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


def fill_graph_bilingual(
    db_path: Path,
    pt_en: dict[str, str],
    pt_zh: dict[str, str],
    soc_en: dict[str, str],
    soc_zh: dict[str, str],
) -> None:
    conn = sqlite3.connect(str(db_path))
    try:
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA synchronous=NORMAL")

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

        node_pt_rows = [
            (pt_en.get(code, name), pt_zh.get(code, ""), code)
            for code, name in conn.execute("SELECT pt_code, pt_name FROM node_pt")
        ]
        conn.executemany(
            "UPDATE node_pt SET pt_name_en=?, pt_name_zh=? WHERE pt_code=?",
            node_pt_rows,
        )

        for table in ["edge_compound_pt_ps", "edge_compound_pt_ss"]:
            edge_rows = [
                (pt_en.get(code, name), pt_zh.get(code, ""), code)
                for code, name in conn.execute(f"SELECT DISTINCT pt_code, pt_name FROM {table}")
            ]
            conn.executemany(
                f"UPDATE {table} SET pt_name_en=?, pt_name_zh=? WHERE pt_code=?",
                edge_rows,
            )

        node_soc_rows = [
            (soc_en.get(code, name), soc_zh.get(code, ""), code)
            for name, code in conn.execute("SELECT soc_name, soc_code FROM node_soc")
        ]
        conn.executemany(
            "UPDATE node_soc SET soc_name_en=?, soc_name_zh=? WHERE soc_code=?",
            node_soc_rows,
        )

        edge_soc_rows = [
            (soc_en.get(code, name), soc_zh.get(code, ""), code)
            for code, name in conn.execute(
                "SELECT DISTINCT soc_code_primary, soc_name_primary FROM edge_pt_soc_primary"
            )
        ]
        conn.executemany(
            "UPDATE edge_pt_soc_primary SET soc_name_en=?, soc_name_zh=? WHERE soc_code_primary=?",
            edge_soc_rows,
        )
        conn.commit()
    finally:
        conn.close()


def build_reference_excels(
    pt_en: dict[str, str],
    pt_zh: dict[str, str],
    soc_en: dict[str, str],
    soc_zh: dict[str, str],
) -> tuple[Path, Path]:
    base = sqlite3.connect(str(BASE_DB))
    try:
        compounds = pd.read_sql_query(
            """
            SELECT
                inchikey,
                smiles_std,
                inchi_std,
                term_count,
                dbid_count,
                dbid_list,
                drugbank_name_list
            FROM node_compound
            ORDER BY inchikey
            """,
            base,
        )
        reactions = pd.read_sql_query(
            """
            SELECT DISTINCT
                p.pt_code,
                p.pt_name AS pt_name_base,
                e.soc_code_primary,
                e.soc_name_primary
            FROM node_pt p
            LEFT JOIN edge_pt_soc_primary e
                ON p.pt_code = e.pt_code
            ORDER BY p.pt_code
            """,
            base,
        )
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

    reactions["pt_name_en"] = reactions.apply(
        lambda r: pt_en.get(r["pt_code"], r["pt_name_base"] if pd.notna(r["pt_name_base"]) else ""),
        axis=1,
    )
    reactions["pt_name_zh"] = reactions["pt_code"].map(lambda x: pt_zh.get(x, ""))
    reactions["soc_name_en"] = reactions["soc_code_primary"].map(lambda x: soc_en.get(x, "") if pd.notna(x) else "")
    reactions["soc_name_zh"] = reactions["soc_code_primary"].map(lambda x: soc_zh.get(x, "") if pd.notna(x) else "")
    for name in ["ge10", "ge20", "ge30"]:
        reactions[f"in_{name}"] = reactions["pt_code"].map(lambda x: 1 if x in graph_pts[name] else 0)

    reactions = reactions[
        [
            "pt_code",
            "pt_name_en",
            "pt_name_zh",
            "soc_code_primary",
            "soc_name_en",
            "soc_name_zh",
            "in_ge10",
            "in_ge20",
            "in_ge30",
        ]
    ]

    compound_xlsx = PACKAGE_DIR / "all_compounds_basic_info.xlsx"
    reaction_xlsx = PACKAGE_DIR / "all_reactions_basic_info_bilingual.xlsx"
    with pd.ExcelWriter(compound_xlsx, engine="openpyxl") as writer:
        compounds.to_excel(writer, sheet_name="compounds", index=False)
    with pd.ExcelWriter(reaction_xlsx, engine="openpyxl") as writer:
        reactions.to_excel(writer, sheet_name="reactions", index=False)
    return compound_xlsx, reaction_xlsx


def sync_graphs_to_package() -> None:
    dst_root = PACKAGE_DIR / "03_pruned_graphs"
    dst_root.mkdir(parents=True, exist_ok=True)
    for name in ["compound_pt_role_core_ge10", "compound_pt_role_core_ge20", "compound_pt_role_core_ge30"]:
        src = STEP14_DIR / name
        dst = dst_root / name
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(src, dst)


def patch_docs() -> None:
    readme = PACKAGE_DIR / "README.md"
    text = readme.read_text(encoding="utf-8")
    addition = """

## 双语反应与基础信息表

1. `all_compounds_basic_info.xlsx`：整理后的全部结构化成分基础信息，包含 `InChIKey`、标准化 `SMILES`、标准化 `InChI`、DrugBank 名称汇总，以及该成分是否进入 `ge10/ge20/ge30` 三个裁剪图。
2. `all_reactions_basic_info_bilingual.xlsx`：整理后的全部反应基础信息，包含 MedDRA `PT` 编码、英文名、中文名、主 `SOC` 编码、主 `SOC` 中英文名称，以及该反应是否进入 `ge10/ge20/ge30` 三个裁剪图。
3. `03_pruned_graphs` 下三个图中的反应和系统节点，现已统一补充中英文双语字段，便于专家快速核对和后续中文写作。
""".strip()
    if "all_reactions_basic_info_bilingual.xlsx" not in text:
        text = text.rstrip() + "\n\n" + addition + "\n"
        readme.write_text(text, encoding="utf-8")

    methods = PACKAGE_DIR / "cleaning_methods_results.md"
    text2 = methods.read_text(encoding="utf-8")
    addition2 = """

## 中英文 MedDRA 补充

1. 在三个裁剪版 FAERS 核心图中，为 `PT` 和 `SOC` 节点统一补充了 MedDRA 英文名称和中文名称。
2. 在 `compound -> PT` 边表中，同步补充了 `PT` 的中英文名称，避免后续专家核对时还要再回查字典。
3. 另外导出两份 Excel 参考表，分别覆盖全部结构化成分和全部反应的基础信息，作为图数据库之外的快速查看入口。
""".strip()
    if "中英文 MedDRA 补充" not in text2:
        text2 = text2.rstrip() + "\n\n" + addition2 + "\n"
        methods.write_text(text2, encoding="utf-8")


def main() -> None:
    pt_en = read_asc_map(MEDDRA_EN / "pt.asc")
    pt_zh = read_asc_map(MEDDRA_ZH / "pt.asc")
    soc_en = read_asc_map(MEDDRA_EN / "soc.asc")
    soc_zh = read_asc_map(MEDDRA_ZH / "soc.asc")

    for db in GRAPH_DIRS.values():
        fill_graph_bilingual(db, pt_en, pt_zh, soc_en, soc_zh)

    build_reference_excels(pt_en, pt_zh, soc_en, soc_zh)
    sync_graphs_to_package()
    patch_docs()
    print("done")


if __name__ == "__main__":
    main()
