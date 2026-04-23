#!/usr/bin/env python3
"""Step4: build InChIKey-reaction dataset from Step2+Step3 (without MedDRA dict).

Core logic:
1) Take Step3 drug terms with non-empty InChI.
2) Normalize InChI with RDKit (cleanup, fragment parent, uncharge, canonical tautomer).
3) Use standardized InChIKey as unique drug structure id.
4) Join with Step2 reactions by primaryid to build InChIKey-PT counts.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

from rdkit import Chem
from rdkit.Chem import inchi as rd_inchi
from rdkit.Chem.MolStandardize import rdMolStandardize


@dataclass
class InchiNormRow:
    inchi_raw: str
    status: str
    inchikey: str
    inchi_std: str
    smiles_std: str
    n_frag_before: int
    n_frag_after: int
    note: str


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build InChIKey-reaction dataset")
    p.add_argument(
        "--step2-db",
        type=Path,
        default=Path(r"D:\博士文件\TCMMKG\data\AEMS_FDA不良反应数据\smiles_adr_project\outputs\step2\faers_step2.sqlite"),
    )
    p.add_argument(
        "--step3-db",
        type=Path,
        default=Path(r"D:\博士文件\TCMMKG\data\AEMS_FDA不良反应数据\smiles_adr_project\outputs\step3\faers_step3.sqlite"),
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path(r"D:\博士文件\TCMMKG\data\AEMS_FDA不良反应数据\smiles_adr_project\outputs\step4"),
    )
    p.add_argument("--min-pair-cases", type=int, default=20)
    p.add_argument(
        "--include-any",
        action="store_true",
        help="Also build all-role (ANY) InChIKey-PT aggregation (much heavier).",
    )
    p.add_argument("--rebuild", action="store_true")
    return p.parse_args()


def scalar(conn: sqlite3.Connection, q: str) -> int:
    return int(conn.execute(q).fetchone()[0])


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


def init_db(conn: sqlite3.Connection, rebuild: bool) -> None:
    conn.execute("PRAGMA journal_mode = WAL;")
    conn.execute("PRAGMA synchronous = NORMAL;")
    conn.execute("PRAGMA temp_store = FILE;")
    cur = conn.cursor()
    if rebuild:
        cur.executescript(
            """
            DROP TABLE IF EXISTS term_inchi_raw;
            DROP TABLE IF EXISTS inchi_norm_map;
            DROP TABLE IF EXISTS term_inchikey_map;
            DROP TABLE IF EXISTS inchikey_meta;
            DROP TABLE IF EXISTS reaction_base;
            DROP TABLE IF EXISTS case_inchikey_any;
            DROP TABLE IF EXISTS case_inchikey_psss;
            DROP TABLE IF EXISTS inchikey_case_counts_any;
            DROP TABLE IF EXISTS inchikey_case_counts_psss;
            DROP TABLE IF EXISTS inchikey_pt_counts_any;
            DROP TABLE IF EXISTS inchikey_pt_counts_psss;
            DROP TABLE IF EXISTS inchikey_pt_psss_model;
            DROP TABLE IF EXISTS inchikey_index;
            DROP TABLE IF EXISTS pt_index;
            DROP TABLE IF EXISTS inchikey_pt_matrix_psss;
            """
        )
    conn.commit()


def normalize_inchi(raw_inchi: str, uncharger: rdMolStandardize.Uncharger, taut_enum: rdMolStandardize.TautomerEnumerator) -> InchiNormRow:
    raw = (raw_inchi or "").strip()
    if not raw:
        return InchiNormRow(raw, "EMPTY", "", "", "", 0, 0, "empty input")

    try:
        mol = rd_inchi.MolFromInchi(raw, sanitize=True, removeHs=True)
    except Exception as exc:
        return InchiNormRow(raw, "PARSE_FAIL", "", "", "", 0, 0, f"MolFromInchi error: {exc}")

    if mol is None:
        return InchiNormRow(raw, "PARSE_FAIL", "", "", "", 0, 0, "MolFromInchi returned None")

    try:
        n_frag_before = len(Chem.GetMolFrags(mol))
        mol = rdMolStandardize.Cleanup(mol)
        mol = rdMolStandardize.FragmentParent(mol)
        mol = uncharger.uncharge(mol)
        mol = rdMolStandardize.Cleanup(mol)
        mol = taut_enum.Canonicalize(mol)
        if mol is None or mol.GetNumAtoms() == 0:
            return InchiNormRow(raw, "EMPTY_PARENT", "", "", "", n_frag_before, 0, "empty parent after standardize")

        n_frag_after = len(Chem.GetMolFrags(mol))
        smiles_std = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        inchi_std = rd_inchi.MolToInchi(mol)
        if not inchi_std:
            return InchiNormRow(raw, "INCHI_FAIL", "", "", smiles_std, n_frag_before, n_frag_after, "MolToInchi empty")
        inchikey = rd_inchi.InchiToInchiKey(inchi_std)
        if not inchikey:
            return InchiNormRow(raw, "INCHIKEY_FAIL", "", inchi_std, smiles_std, n_frag_before, n_frag_after, "InchiToInchiKey empty")

        note = "ok"
        if n_frag_before > 1:
            note = f"fragment_parent_from_{n_frag_before}"
        return InchiNormRow(raw, "OK", inchikey, inchi_std, smiles_std, n_frag_before, n_frag_after, note)
    except Exception as exc:
        return InchiNormRow(raw, "STANDARDIZE_FAIL", "", "", "", 0, 0, f"standardize error: {exc}")


def build_inchi_norm_map(conn: sqlite3.Connection) -> Tuple[int, int]:
    cur = conn.cursor()
    raw_inchis = [r[0] for r in cur.execute("SELECT DISTINCT inchi_raw FROM term_inchi_raw WHERE TRIM(COALESCE(inchi_raw,''))<>''")]

    uncharger = rdMolStandardize.Uncharger()
    taut_enum = rdMolStandardize.TautomerEnumerator()

    norm_rows: List[Tuple[str, str, str, str, str, int, int, str]] = []
    ok = 0
    for raw in raw_inchis:
        nr = normalize_inchi(raw, uncharger, taut_enum)
        if nr.status == "OK":
            ok += 1
        norm_rows.append(
            (
                nr.inchi_raw,
                nr.status,
                nr.inchikey,
                nr.inchi_std,
                nr.smiles_std,
                nr.n_frag_before,
                nr.n_frag_after,
                nr.note,
            )
        )

    cur.execute("DROP TABLE IF EXISTS inchi_norm_map")
    cur.execute(
        """
        CREATE TABLE inchi_norm_map (
            inchi_raw TEXT PRIMARY KEY,
            status TEXT NOT NULL,
            inchikey TEXT,
            inchi_std TEXT,
            smiles_std TEXT,
            n_frag_before INTEGER,
            n_frag_after INTEGER,
            note TEXT
        )
        """
    )
    cur.executemany(
        """
        INSERT INTO inchi_norm_map(
            inchi_raw,status,inchikey,inchi_std,smiles_std,n_frag_before,n_frag_after,note
        ) VALUES (?,?,?,?,?,?,?,?)
        """,
        norm_rows,
    )
    cur.executescript(
        """
        CREATE INDEX idx_inchi_norm_status ON inchi_norm_map(status);
        CREATE INDEX idx_inchi_norm_inchikey ON inchi_norm_map(inchikey);
        """
    )
    conn.commit()
    return len(raw_inchis), ok


def main() -> int:
    args = parse_args()
    started = dt.datetime.now().isoformat(timespec="seconds")

    if not args.step2_db.exists():
        print(f"[ERROR] step2 db not found: {args.step2_db}")
        return 1
    if not args.step3_db.exists():
        print(f"[ERROR] step3 db not found: {args.step3_db}")
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    step4_db = args.output_dir / "faers_step4.sqlite"

    conn = sqlite3.connect(step4_db)
    try:
        init_db(conn, rebuild=args.rebuild)
        conn.execute("ATTACH DATABASE ? AS s2", (str(args.step2_db),))
        conn.execute("ATTACH DATABASE ? AS s3", (str(args.step3_db),))

        print("[INFO] extracting term_inchi_raw from step3 ...", flush=True)
        conn.executescript(
            """
            DROP TABLE IF EXISTS term_inchi_raw;
            CREATE TABLE term_inchi_raw AS
            SELECT
                term_std,
                term_norm,
                n_rows,
                n_cases,
                dbid,
                drugbank_name,
                inchi AS inchi_raw,
                mapping_stage,
                manual_action,
                manual_target
            FROM s3.drug_term_final
            WHERE TRIM(COALESCE(inchi,'')) <> ''
              AND UPPER(TRIM(COALESCE(mapping_stage,''))) <> 'MANUAL_EXCLUDE';

            CREATE INDEX idx_term_inchi_raw_term ON term_inchi_raw(term_std);
            CREATE INDEX idx_term_inchi_raw_inchi ON term_inchi_raw(inchi_raw);
            """
        )
        conn.commit()

        print("[INFO] RDKit normalizing InChI -> InChIKey/SMILES ...", flush=True)
        raw_inchi_total, norm_ok_total = build_inchi_norm_map(conn)

        print("[INFO] building term_inchikey_map and inchikey_meta ...", flush=True)
        conn.executescript(
            """
            DROP TABLE IF EXISTS term_inchikey_map;
            CREATE TABLE term_inchikey_map AS
            SELECT
                t.term_std,
                t.term_norm,
                t.n_rows,
                t.n_cases,
                t.dbid,
                t.drugbank_name,
                t.inchi_raw,
                n.inchikey,
                n.inchi_std,
                n.smiles_std,
                n.n_frag_before,
                n.n_frag_after,
                n.note AS norm_note,
                t.mapping_stage,
                t.manual_action,
                t.manual_target
            FROM term_inchi_raw t
            INNER JOIN inchi_norm_map n
                ON t.inchi_raw = n.inchi_raw
            WHERE n.status = 'OK'
              AND TRIM(COALESCE(n.inchikey,'')) <> '';

            CREATE INDEX idx_term_inchikey_term ON term_inchikey_map(term_std);
            CREATE INDEX idx_term_inchikey_key ON term_inchikey_map(inchikey);

            DROP TABLE IF EXISTS inchikey_meta;
            CREATE TABLE inchikey_meta AS
            SELECT
                inchikey,
                MIN(inchi_std) AS inchi_std,
                MIN(smiles_std) AS smiles_std,
                COUNT(DISTINCT term_std) AS term_count,
                COUNT(DISTINCT dbid) AS dbid_count,
                GROUP_CONCAT(DISTINCT dbid) AS dbid_list,
                GROUP_CONCAT(DISTINCT drugbank_name) AS drugbank_name_list
            FROM term_inchikey_map
            GROUP BY inchikey;

            CREATE INDEX idx_inchikey_meta_key ON inchikey_meta(inchikey);
            """
        )
        conn.commit()

        print("[INFO] preparing reaction base ...", flush=True)
        conn.executescript(
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
            """
        )
        conn.commit()

        print("[INFO] building case_inchikey_any/psss ...", flush=True)
        conn.executescript(
            """
            DROP TABLE IF EXISTS case_inchikey_any;
            CREATE TABLE case_inchikey_any AS
            SELECT DISTINCT
                d.primaryid,
                d.caseid,
                m.inchikey
            FROM s2.case_drug d
            INNER JOIN term_inchikey_map m
                ON d.drug_key_std = m.term_std
            WHERE TRIM(COALESCE(d.primaryid,'')) <> ''
              AND TRIM(COALESCE(m.inchikey,'')) <> '';

            CREATE INDEX idx_case_inchikey_any_pid ON case_inchikey_any(primaryid);
            CREATE INDEX idx_case_inchikey_any_key ON case_inchikey_any(inchikey);

            DROP TABLE IF EXISTS case_inchikey_psss;
            CREATE TABLE case_inchikey_psss AS
            SELECT DISTINCT
                d.primaryid,
                d.caseid,
                m.inchikey
            FROM s2.case_drug d
            INNER JOIN term_inchikey_map m
                ON d.drug_key_std = m.term_std
            WHERE TRIM(COALESCE(d.primaryid,'')) <> ''
              AND TRIM(COALESCE(m.inchikey,'')) <> ''
              AND UPPER(TRIM(COALESCE(d.role_cod,''))) IN ('PS','SS');

            CREATE INDEX idx_case_inchikey_psss_pid ON case_inchikey_psss(primaryid);
            CREATE INDEX idx_case_inchikey_psss_key ON case_inchikey_psss(inchikey);
            """
        )
        conn.commit()

        print("[INFO] aggregating inchikey counts ...", flush=True)
        conn.executescript(
            """
            DROP TABLE IF EXISTS inchikey_case_counts_any;
            CREATE TABLE inchikey_case_counts_any AS
            SELECT inchikey, COUNT(*) AS n_cases
            FROM case_inchikey_any
            GROUP BY inchikey;
            CREATE INDEX idx_ik_case_any ON inchikey_case_counts_any(inchikey);

            DROP TABLE IF EXISTS inchikey_case_counts_psss;
            CREATE TABLE inchikey_case_counts_psss AS
            SELECT inchikey, COUNT(*) AS n_cases
            FROM case_inchikey_psss
            GROUP BY inchikey;
            CREATE INDEX idx_ik_case_psss ON inchikey_case_counts_psss(inchikey);
            """
        )
        conn.commit()

        if args.include_any:
            print("[INFO] aggregating inchikey-PT counts (ANY) ...", flush=True)
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

        print("[INFO] aggregating inchikey-PT counts (PS/SS) ...", flush=True)
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

        print("[INFO] building model tables/matrix ...", flush=True)
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

        print("[INFO] exporting step4 files ...", flush=True)
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
            "min_pair_cases": args.min_pair_cases,
            "include_any": bool(args.include_any),
            "term_inchi_raw_rows": scalar(conn, "SELECT COUNT(*) FROM term_inchi_raw"),
            "term_inchi_raw_distinct_terms": scalar(conn, "SELECT COUNT(DISTINCT term_std) FROM term_inchi_raw"),
            "term_inchi_raw_distinct_inchi": scalar(conn, "SELECT COUNT(DISTINCT inchi_raw) FROM term_inchi_raw"),
            "raw_inchi_total": raw_inchi_total,
            "raw_inchi_norm_ok": norm_ok_total,
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
        print("[DONE] Step4 completed", flush=True)
        print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    raise SystemExit(main())
