"""Step 6d: Build pruned core graph releases.

Purpose:
- derive lighter public graph releases from the full Step6 core graph
- keep only compound-PT edges above fixed `n_reports` thresholds
- export ge10, ge20, and ge30 graph packages with summary files

Inputs:
- outputs/step6_core_graphs/core_graph_full/step6_compound_pt_soc_core_full.sqlite

Outputs:
- outputs/step6_core_graphs/core_graph_ge10/
- outputs/step6_core_graphs/core_graph_ge20/
- outputs/step6_core_graphs/core_graph_ge30/
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC_DB = ROOT / "outputs" / "step6_core_graphs" / "core_graph_full" / "step6_compound_pt_soc_core_full.sqlite"
BASE_OUT = ROOT / "outputs" / "step6_core_graphs"
THRESHOLDS = [10, 20, 30]


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def scalar(conn: sqlite3.Connection, q: str) -> int:
    return int(conn.execute(q).fetchone()[0])


def export_csv(conn: sqlite3.Connection, q: str, out_path: Path) -> None:
    cur = conn.cursor()
    cur.execute(q)
    cols = [d[0] for d in cur.description]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        f.write(",".join(cols) + "\n")
        for row in cur:
            vals = []
            for v in row:
                s = "" if v is None else str(v)
                if "," in s or "\"" in s:
                    s = '"' + s.replace('"', '""') + '"'
                vals.append(s)
            f.write(",".join(vals) + "\n")


def build_one(threshold: int) -> dict[str, object]:
    out_dir = BASE_OUT / f"step6_core_graph_ge{threshold}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_db = out_dir / f"step6_compound_pt_soc_core_ge{threshold}.sqlite"
    if out_db.exists():
        out_db.unlink()

    src = sqlite3.connect(str(SRC_DB))
    dst = sqlite3.connect(str(out_db))
    dst.execute("PRAGMA journal_mode=DELETE;")
    dst.execute("PRAGMA synchronous=OFF;")
    src_attach = str(SRC_DB).replace("\\", "/")
    dst.execute(f"ATTACH DATABASE '{src_attach}' AS src;")

    try:
        dst.executescript(
            """
            CREATE TABLE edge_compound_pt_ps (
                inchikey TEXT NOT NULL,
                pt_code TEXT NOT NULL,
                pt_name TEXT,
                n_reports INTEGER NOT NULL,
                PRIMARY KEY (inchikey, pt_code)
            );
            CREATE TABLE edge_compound_pt_ss (
                inchikey TEXT NOT NULL,
                pt_code TEXT NOT NULL,
                pt_name TEXT,
                n_reports INTEGER NOT NULL,
                PRIMARY KEY (inchikey, pt_code)
            );
            CREATE TABLE node_compound (
                inchikey TEXT,
                inchi_std TEXT,
                smiles_std TEXT,
                term_count INTEGER,
                dbid_count INTEGER,
                dbid_list TEXT,
                drugbank_name_list TEXT
            );
            CREATE TABLE node_pt (
                pt_code TEXT,
                pt_name TEXT
            );
            CREATE TABLE node_soc (
                soc_name TEXT,
                soc_code TEXT
            );
            CREATE TABLE edge_pt_soc_primary (
                pt_code TEXT,
                soc_code_primary TEXT,
                soc_name_primary TEXT
            );
            """
        )

        for role in ["ps", "ss"]:
            dst.execute(
                f"""
                INSERT INTO main.edge_compound_pt_{role}
                SELECT inchikey, pt_code, pt_name, n_reports
                FROM src.edge_compound_pt_{role}
                WHERE n_reports >= ?
                """,
                (threshold,),
            )

        dst.execute(
            """
            INSERT INTO node_compound
            SELECT *
            FROM src.node_compound
            WHERE inchikey IN (
                SELECT DISTINCT inchikey FROM edge_compound_pt_ps
                UNION
                SELECT DISTINCT inchikey FROM edge_compound_pt_ss
            )
            """
        )
        dst.execute(
            """
            INSERT INTO node_pt
            SELECT *
            FROM src.node_pt
            WHERE pt_code IN (
                SELECT DISTINCT pt_code FROM edge_compound_pt_ps
                UNION
                SELECT DISTINCT pt_code FROM edge_compound_pt_ss
            )
            """
        )
        dst.execute(
            """
            INSERT INTO edge_pt_soc_primary
            SELECT *
            FROM src.edge_pt_soc_primary
            WHERE pt_code IN (SELECT pt_code FROM node_pt)
            """
        )
        dst.execute(
            """
            INSERT INTO node_soc
            SELECT *
            FROM src.node_soc
            WHERE soc_code IN (
                SELECT DISTINCT soc_code_primary FROM edge_pt_soc_primary
            )
            """
        )
        dst.commit()
        dst.execute("VACUUM")
        dst.commit()

        node_summary = {
            "compound": scalar(dst, "SELECT COUNT(*) FROM node_compound"),
            "pt": scalar(dst, "SELECT COUNT(*) FROM node_pt"),
            "soc": scalar(dst, "SELECT COUNT(*) FROM node_soc"),
        }
        edge_summary = {
            "compound_has_pt_ps": scalar(dst, "SELECT COUNT(*) FROM edge_compound_pt_ps"),
            "compound_has_pt_ss": scalar(dst, "SELECT COUNT(*) FROM edge_compound_pt_ss"),
            "pt_belongs_to_primary_soc": scalar(dst, "SELECT COUNT(*) FROM edge_pt_soc_primary"),
        }
        weight_summary = {
            "ps_total_reports": scalar(dst, "SELECT COALESCE(SUM(n_reports),0) FROM edge_compound_pt_ps"),
            "ss_total_reports": scalar(dst, "SELECT COALESCE(SUM(n_reports),0) FROM edge_compound_pt_ss"),
        }
    finally:
        src.close()
        dst.close()

    size_bytes = out_db.stat().st_size
    summary = {
        "threshold": threshold,
        "sqlite": str(out_db),
        "size_bytes": size_bytes,
        "size_mb": round(size_bytes / (1024 * 1024), 2),
        "nodes": node_summary,
        "edges": edge_summary,
        "weights": weight_summary,
    }

    write_text(out_dir / "graph_summary.json", json.dumps(summary, ensure_ascii=False, indent=2))
    write_text(
        out_dir / "README.md",
        "\n".join(
            [
                f"# FAERS Core Graph (n_reports >= {threshold})",
                "",
                "这是在 `core_graph_full` 基础上裁剪得到的简化版主图。",
                "",
                "保留规则：",
                f"1. `edge_compound_pt_ps`: `n_reports >= {threshold}`",
                f"2. `edge_compound_pt_ss`: `n_reports >= {threshold}`",
                "3. `edge_pt_soc_primary`: 仅保留仍然被保留 PT 连接到的层级边",
                "",
                "节点：`compound`、`pt`、`soc`",
            ]
        ),
    )

    # Export small summaries
    conn = sqlite3.connect(str(out_db))
    try:
        export_csv(
            conn,
            f"""
            SELECT 'compound' AS node_type, COUNT(*) AS n_nodes FROM node_compound
            UNION ALL
            SELECT 'pt', COUNT(*) FROM node_pt
            UNION ALL
            SELECT 'soc', COUNT(*) FROM node_soc
            """,
            out_dir / "node_summary.csv",
        )
        export_csv(
            conn,
            f"""
            SELECT 'compound_has_pt_ps' AS edge_type, 'compound' AS src_type, 'pt' AS dst_type, COUNT(*) AS n_edges FROM edge_compound_pt_ps
            UNION ALL
            SELECT 'compound_has_pt_ss', 'compound', 'pt', COUNT(*) FROM edge_compound_pt_ss
            UNION ALL
            SELECT 'pt_belongs_to_primary_soc', 'pt', 'soc', COUNT(*) FROM edge_pt_soc_primary
            """,
            out_dir / "edge_summary.csv",
        )
    finally:
        conn.close()

    return summary


def main() -> None:
    summaries = []
    for threshold in THRESHOLDS:
        summaries.append(build_one(threshold))
    print(json.dumps(summaries, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
