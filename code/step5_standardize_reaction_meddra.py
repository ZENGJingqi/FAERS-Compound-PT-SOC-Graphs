#!/usr/bin/env python3
"""Step5: MedDRA standardization for reactions using existing Step2 database.

This script does an incremental standardization pass:
- Read distinct `pt_norm` from Step2 `case_reaction`.
- Map to MedDRA v29 PT code/name by PT exact, LLT->PT, and curated manual code hints.
- Build Step5 sqlite outputs with mapping tables and mapped reaction base for downstream Step4 rebuild.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import re
import sqlite3
import unicodedata
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple


# Reuse the curated manual hints used in Step2.
REACTION_MANUAL_CODE_HINTS = {
    "ANO-RECTAL STENOSIS": "10002581",
    "HER-2 POSITIVE BREAST CANCER": "10065430",
    "STAPHYLOCOCCAL IDENTIFICATION TEST POSITIVE": "10067140",
    "STREPTOCOCCAL IDENTIFICATION TEST": "10067006",
    "STREPTOCOCCAL SEROLOGY": "10059987",
    "ANO-RECTAL ULCER": "10002582",
    "BLASTIC PLASMACYTOID DENDRITRIC CELL NEOPLASIA": "10075460",
    "CORNELIA DE-LANGE SYNDROME": "10077707",
    "FRONTAL SINUS OPERATIONS": "10017379",
    "HER-2 POSITIVE GASTRIC CANCER": "10066896",
    "MAXILLARY ANTRUM OPERATIONS": "10026950",
    "METHICILLIN-RESISTANT STAPHYLOCOCCAL AUREUS TEST NEGATIVE": "10053428",
    "METHICILLIN-RESISTANT STAPHYLOCOCCAL AUREUS TEST POSITIVE": "10053427",
    "PAROVARIAN CYST": "10052456",
    "STAPHYLOCOCCAL IDENTIFICATION TEST NEGATIVE": "10067005",
    "STREPTOCOCCAL IDENTIFICATION TEST POSITIVE": "10067004",
    "AEROMONA INFECTION": "10054205",
    "DISBACTERIOSIS": "10064389",
    "GASTRO-ENTEROSTOMY": "10017873",
    "HER-2 PROTEIN OVEREXPRESSION": "10075638",
    "HYPOTHALAMO-PITUITARY DISORDERS": "10021111",
    "SUPERIOR VENA CAVAL OCCLUSION": "10058988",
    "CAPNOCYTOPHAGIA INFECTION": "10061738",
    "EAGLES SYNDROME": "10066835",
    "EVAN'S SYNDROME": "10053873",
    "GASTRO-INTESTINAL FISTULA": "10071258",
    "GLYCOPEPTIDE ANTIBIOTIC RESISTANT STAPHYLOCOCCAL AUREUS INFECTION": "10052101",
    "HEPATOBILLIARY DISORDER PROPHYLAXIS": "10081385",
    "MENINGEOMAS SURGERY": "10053765",
    "METHICILLIN-RESISTANT STAPHYLOCOCCAL AUREUS TEST": "10053429",
    "SPHENOID SINUS OPERATIONS": "10041508",
    "STREPTOCOCCAL SEROLOGY POSITIVE": "10059988",
    "SUPERIOR VENA CAVAL STENOSIS": "10064771",
    "GASTRO-JEJUNOSTOMY": "10017882",
    "IMMUNE-MEDIATED ADRENAL INSUFICIENCY": "10085547",
    "ORAL APPLIANCE": "10085270",
    "VAGINAL RING": "10082353",
    "VAGINAL CUFF": "10088846",
    "MOREL-LAVELLEE SEROMA": "10088873",
    "RABSON MENDENHALL SYNDROME": "10088742",
}


@dataclass
class MedDRAResources:
    pt_name_to_code: Dict[str, str]
    pt_code_to_name: Dict[str, str]
    llt_name_to_pt_code: Dict[str, str]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Step5 MedDRA reaction standardization")
    p.add_argument("--step2-db", type=Path, default=Path("outputs/step2/faers_step2.sqlite"))
    p.add_argument(
        "--meddra-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "reference_data" / "MedDRA" / "MedDRA_29_0_English",
        help="Path containing MedAscii/pt.asc and MedAscii/llt.asc",
    )
    p.add_argument("--output-dir", type=Path, default=Path("outputs/step5"))
    p.add_argument("--rebuild", action="store_true")
    return p.parse_args()


def normalize_text(s: Optional[str], upper: bool = True) -> str:
    if s is None:
        return ""
    t = unicodedata.normalize("NFKC", s)
    t = t.replace("\ufeff", "")
    t = t.replace("\u2013", "-").replace("\u2014", "-")
    t = t.replace("\u2018", "'").replace("\u2019", "'")
    t = t.replace("\u201c", '"').replace("\u201d", '"')
    t = t.strip()
    t = re.sub(r"\s+", " ", t)
    if upper:
        t = t.upper()
    return t


def normalize_reaction_term(s: Optional[str]) -> str:
    t = normalize_text(s, upper=True)
    if not t:
        return ""
    t = re.sub(r"\s*[-]+\s*", "-", t)
    t = re.sub(r"\s*/\s*", " / ", t)
    t = re.sub(r"\s+", " ", t).strip()
    return t


def find_medascii_dir(meddra_dir: Path) -> Optional[Path]:
    if not meddra_dir.exists():
        return None
    direct = meddra_dir / "MedAscii"
    if direct.exists() and direct.is_dir():
        return direct
    for d in meddra_dir.rglob("MedAscii"):
        if d.is_dir():
            return d
    return None


def load_meddra_resources(meddra_dir: Path) -> Tuple[MedDRAResources, str]:
    medascii = find_medascii_dir(meddra_dir)
    if medascii is None:
        raise FileNotFoundError(f"MedAscii not found under: {meddra_dir}")

    pt_path = medascii / "pt.asc"
    llt_path = medascii / "llt.asc"
    rel_path = medascii / "meddra_release.asc"
    if not pt_path.exists() or not llt_path.exists():
        raise FileNotFoundError(f"pt.asc/llt.asc missing under: {medascii}")

    pt_name_to_code: Dict[str, str] = {}
    pt_code_to_name: Dict[str, str] = {}
    llt_name_to_pt_code: Dict[str, str] = {}

    with pt_path.open("r", encoding="latin-1", errors="replace", newline="") as f:
        for line in f:
            row = line.rstrip("\r\n").split("$")
            if len(row) < 2:
                continue
            pt_code = row[0].strip()
            pt_name = row[1].strip()
            if not pt_code or not pt_name:
                continue
            k = normalize_reaction_term(pt_name)
            if k and k not in pt_name_to_code:
                pt_name_to_code[k] = pt_code
            if pt_code not in pt_code_to_name:
                pt_code_to_name[pt_code] = pt_name

    with llt_path.open("r", encoding="latin-1", errors="replace", newline="") as f:
        for line in f:
            row = line.rstrip("\r\n").split("$")
            if len(row) < 3:
                continue
            llt_name = row[1].strip()
            pt_code = row[2].strip()
            if not llt_name or not pt_code:
                continue
            k = normalize_reaction_term(llt_name)
            if k and k not in llt_name_to_pt_code:
                llt_name_to_pt_code[k] = pt_code

    release = ""
    if rel_path.exists():
        release = rel_path.read_text(encoding="utf-8", errors="ignore").strip()

    return MedDRAResources(pt_name_to_code, pt_code_to_name, llt_name_to_pt_code), release


def map_reaction_term(pt_norm: str, med: MedDRAResources) -> Tuple[str, str, str, str]:
    t = normalize_reaction_term(pt_norm)
    hint_code = REACTION_MANUAL_CODE_HINTS.get(t, "")

    if not t:
        return t, hint_code, "", "", "EMPTY_PT"

    pt_code = med.pt_name_to_code.get(t, "")
    if pt_code:
        return t, hint_code, pt_code, med.pt_code_to_name.get(pt_code, t), "PT_EXACT"

    llt_to_pt = med.llt_name_to_pt_code.get(t, "")
    if llt_to_pt:
        return t, hint_code, llt_to_pt, med.pt_code_to_name.get(llt_to_pt, ""), "LLT_TO_PT"

    if hint_code and hint_code in med.pt_code_to_name:
        return t, hint_code, hint_code, med.pt_code_to_name[hint_code], "MANUAL_CODE_HINT"

    if hint_code:
        return t, hint_code, "", "", "MANUAL_HINT_CODE_UNRESOLVED"

    return t, "", "", "", "UNRESOLVED_WITH_MEDDRA"


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

    if not args.step2_db.exists():
        print(f"[ERROR] step2 db not found: {args.step2_db}", flush=True)
        return 1

    print("[INFO] loading MedDRA resources ...", flush=True)
    med, release = load_meddra_resources(args.meddra_dir)
    print(
        f"[INFO] MedDRA loaded: PT={len(med.pt_name_to_code)}, LLT={len(med.llt_name_to_pt_code)}, release='{release}'",
        flush=True,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    step5_db = args.output_dir / "faers_step5.sqlite"
    conn = sqlite3.connect(step5_db)

    try:
        conn.execute("PRAGMA journal_mode = WAL;")
        conn.execute("PRAGMA synchronous = NORMAL;")
        conn.execute("PRAGMA temp_store = FILE;")
        conn.execute("ATTACH DATABASE ? AS s2", (str(args.step2_db),))

        if args.rebuild:
            conn.executescript(
                """
                DROP TABLE IF EXISTS reaction_term_meddra_map;
                DROP TABLE IF EXISTS reaction_unresolved_top;
                DROP TABLE IF EXISTS reaction_base_meddra;
                DROP TABLE IF EXISTS pt_case_counts_meddra;
                """
            )
            conn.commit()

        print("[INFO] collecting distinct pt_norm from step2.case_reaction ...", flush=True)
        terms = [
            r[0]
            for r in conn.execute(
                "SELECT DISTINCT pt_norm FROM s2.case_reaction WHERE TRIM(COALESCE(pt_norm,'')) <> ''"
            )
        ]
        print(f"[INFO] distinct reaction terms: {len(terms)}", flush=True)

        print("[INFO] building reaction_term_meddra_map ...", flush=True)
        mapped_rows = []
        source_counter = Counter()
        for term in terms:
            t, hint, code, name, source = map_reaction_term(term, med)
            source_counter[source] += 1
            mapped_rows.append((t, hint, code, name, source))

        conn.executescript(
            """
            DROP TABLE IF EXISTS reaction_term_meddra_map;
            CREATE TABLE reaction_term_meddra_map (
                pt_norm TEXT PRIMARY KEY,
                pt_manual_code_hint TEXT,
                pt_meddra_code TEXT,
                pt_meddra_pt TEXT,
                pt_mapping_source TEXT
            );
            """
        )
        conn.executemany(
            """
            INSERT INTO reaction_term_meddra_map(
                pt_norm,pt_manual_code_hint,pt_meddra_code,pt_meddra_pt,pt_mapping_source
            ) VALUES (?,?,?,?,?)
            """,
            mapped_rows,
        )
        conn.executescript(
            """
            CREATE INDEX idx_map_source ON reaction_term_meddra_map(pt_mapping_source);
            CREATE INDEX idx_map_code ON reaction_term_meddra_map(pt_meddra_code);
            """
        )
        conn.commit()

        print("[INFO] building reaction_unresolved_top ...", flush=True)
        conn.executescript(
            """
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
            """
        )
        conn.commit()

        print("[INFO] building reaction_base_meddra (mapped only) ...", flush=True)
        conn.executescript(
            """
            DROP TABLE IF EXISTS reaction_base_meddra;
            CREATE TABLE reaction_base_meddra AS
            SELECT DISTINCT
                c.primaryid,
                c.caseid,
                m.pt_meddra_code AS pt_code,
                m.pt_meddra_pt AS pt_name,
                m.pt_mapping_source AS map_source
            FROM s2.case_reaction c
            INNER JOIN reaction_term_meddra_map m
                ON c.pt_norm = m.pt_norm
            WHERE TRIM(COALESCE(m.pt_meddra_code,'')) <> '';

            CREATE INDEX idx_reaction_base_meddra_pid ON reaction_base_meddra(primaryid);
            CREATE INDEX idx_reaction_base_meddra_ptcode ON reaction_base_meddra(pt_code);
            CREATE INDEX idx_reaction_base_meddra_ptname ON reaction_base_meddra(pt_name);
            """
        )
        conn.commit()

        print("[INFO] building pt_case_counts_meddra ...", flush=True)
        conn.executescript(
            """
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

        print("[INFO] exporting step5 artifacts ...", flush=True)
        n_map = export_csv(
            conn,
            "SELECT * FROM reaction_term_meddra_map ORDER BY pt_mapping_source, pt_norm",
            args.output_dir / "reaction_term_meddra_map.csv",
        )
        n_unres = export_csv(
            conn,
            "SELECT * FROM reaction_unresolved_top ORDER BY cnt DESC, pt_norm",
            args.output_dir / "reaction_unresolved_top.csv",
        )
        n_unres_top = export_csv(
            conn,
            "SELECT * FROM reaction_unresolved_top ORDER BY cnt DESC, pt_norm LIMIT 1000",
            args.output_dir / "reaction_unresolved_top1000.csv",
        )
        n_pt_counts = export_csv(
            conn,
            "SELECT * FROM pt_case_counts_meddra ORDER BY n_primaryid DESC, pt_code",
            args.output_dir / "pt_case_counts_meddra.csv",
        )

        finished = dt.datetime.now().isoformat(timespec="seconds")
        mapped_terms = scalar(
            conn,
            "SELECT COUNT(*) FROM reaction_term_meddra_map WHERE TRIM(COALESCE(pt_meddra_code,'')) <> ''",
        )
        total_terms = scalar(conn, "SELECT COUNT(*) FROM reaction_term_meddra_map")
        unresolved_terms = total_terms - mapped_terms
        total_rows = scalar(conn, "SELECT COUNT(*) FROM s2.case_reaction")
        mapped_rows_count = scalar(
            conn,
            """
            SELECT COUNT(*)
            FROM s2.case_reaction c
            INNER JOIN reaction_term_meddra_map m ON c.pt_norm = m.pt_norm
            WHERE TRIM(COALESCE(m.pt_meddra_code,'')) <> ''
            """,
        )

        report = {
            "started_at": started,
            "finished_at": finished,
            "meddra_dir": str(args.meddra_dir),
            "meddra_release": release,
            "step2_db": str(args.step2_db),
            "step5_db": str(step5_db),
            "reaction_rows_total": total_rows,
            "reaction_rows_mapped": mapped_rows_count,
            "reaction_rows_unmapped": total_rows - mapped_rows_count,
            "reaction_row_coverage_mapped": (mapped_rows_count / total_rows) if total_rows else 0.0,
            "term_total": total_terms,
            "term_mapped": mapped_terms,
            "term_unmapped": unresolved_terms,
            "term_coverage_mapped": (mapped_terms / total_terms) if total_terms else 0.0,
            "pt_distinct_mapped": scalar(conn, "SELECT COUNT(DISTINCT pt_code) FROM reaction_base_meddra"),
            "primaryid_distinct_mapped": scalar(
                conn, "SELECT COUNT(DISTINCT primaryid) FROM reaction_base_meddra"
            ),
            "mapping_source_counts": dict(source_counter),
            "exports": {
                "reaction_term_meddra_map_rows": n_map,
                "reaction_unresolved_rows": n_unres,
                "reaction_unresolved_top1000_rows": n_unres_top,
                "pt_case_counts_rows": n_pt_counts,
            },
        }

        (args.output_dir / "step5_report.json").write_text(
            json.dumps(report, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

        conn.execute("PRAGMA wal_checkpoint(TRUNCATE);")
        print("[DONE] Step5 MedDRA standardization completed", flush=True)
        print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    raise SystemExit(main())
