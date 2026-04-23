#!/usr/bin/env python3
"""Step 2: Build standardized case-drug and case-reaction datasets.

Design highlights:
- Filter DRUG/REAC rows using Step1 selected latest case keys.
- Standardize drug/reaction strings conservatively.
- Optional MedDRA PT/LLT mapping if a local MedDRA directory is provided.
- Persist into SQLite with dedup constraints, then export CSVs.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import io
import json
import re
import sqlite3
import sys
import unicodedata
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

ZIP_RE = re.compile(r"(?i)^(?:f?aers)_ascii_(20\d{2})[qQ]([1-4])\.zip$")
DRUG_ENTRY_RE = re.compile(r"(?i)(^|/)drug\d{2}q[1-4](?:_new)?\.txt$")
REAC_ENTRY_RE = re.compile(r"(?i)(^|/)reac\d{2}q[1-4](?:_new)?\.txt$")

KNOWN_SPLIT_FIXES = {
    (2011, 2, "drug"): "7475791$1016572490$SS$DOXORUBICIN",
    (2011, 3, "drug"): "7652730$1017255397$SS$BEVACIZUMAB",
    (2011, 4, "drug"): "7941354$1018188213$SS$MEMANTINE HYDROCHLORIDE",
}

# From faers R package clean_reac_pt() curated mapping hints (term -> MedDRA code)
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

DRUG_SYNONYMS = {
    "primaryid": ("primaryid", "isr"),
    "caseid": ("caseid", "case"),
    "drug_seq": ("drug_seq",),
    "role_cod": ("role_cod",),
    "drugname": ("drugname",),
    "prod_ai": ("prod_ai",),
    "route": ("route",),
    "dose_vbm": ("dose_vbm",),
    "val_vbm": ("val_vbm",),
}

REAC_SYNONYMS = {
    "primaryid": ("primaryid", "isr"),
    "caseid": ("caseid", "case"),
    "pt": ("pt",),
    "drug_rec_act": ("drug_rec_act",),
}


@dataclass(frozen=True)
class ZipInfo:
    path: Path
    year: int
    quarter: int


@dataclass
class FileParseStats:
    raw_rows: int = 0
    filtered_rows: int = 0
    inserted_rows: int = 0
    skipped_no_primaryid: int = 0
    trailing_trimmed: int = 0
    padded_short: int = 0
    merged_long: int = 0
    caseid_overridden: int = 0


@dataclass
class MedDRAResources:
    loaded: bool
    pt_name_to_code: Dict[str, str]
    pt_code_to_name: Dict[str, str]
    llt_name_to_pt_code: Dict[str, str]


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


def normalize_drug_term(s: Optional[str]) -> str:
    t = normalize_text(s, upper=True)
    if not t:
        return ""
    t = t.strip("\"'")
    t = re.sub(r"\s*([,;/+()\[\]-])\s*", r" \1 ", t)
    t = re.sub(r"\s+", " ", t).strip()
    return t


def normalize_reaction_term(s: Optional[str]) -> str:
    t = normalize_text(s, upper=True)
    if not t:
        return ""
    t = re.sub(r"\s*[-]+\s*", "-", t)
    t = re.sub(r"\s*/\s*", " / ", t)
    t = re.sub(r"\s+", " ", t).strip()
    return t


def collect_zip_files(data_root: Path) -> List[ZipInfo]:
    out: List[ZipInfo] = []
    for p in data_root.glob("*aers_ascii_*.zip"):
        m = ZIP_RE.match(p.name)
        if not m:
            continue
        out.append(ZipInfo(path=p, year=int(m.group(1)), quarter=int(m.group(2))))
    out.sort(key=lambda x: (x.year, x.quarter, x.path.name.lower()))
    return out


def find_entry(archive: zipfile.ZipFile, pattern: re.Pattern[str]) -> Optional[zipfile.ZipInfo]:
    for e in archive.infolist():
        if pattern.search(e.filename):
            return e
    return None


def build_header_index(header: Sequence[str], synonyms: Dict[str, Tuple[str, ...]]) -> Dict[str, Optional[int]]:
    lowered = [normalize_text(h, upper=True).lower() for h in header]
    idx: Dict[str, Optional[int]] = {}
    for key, aliases in synonyms.items():
        found = None
        for a in aliases:
            if a in lowered:
                found = lowered.index(a)
                break
        idx[key] = found
    return idx


def cell(row: Sequence[str], idx: Optional[int]) -> str:
    if idx is None:
        return ""
    if idx < 0 or idx >= len(row):
        return ""
    return row[idx].strip()


def iter_fixed_lines(
    raw: io.TextIOWrapper,
    year: int,
    quarter: int,
    table_kind: str,
) -> Iterator[str]:
    split_token = KNOWN_SPLIT_FIXES.get((year, quarter, table_kind))
    for line in raw:
        if split_token and split_token in line and not line.startswith(split_token):
            line = line.replace(split_token, "\n" + split_token)
        for sub in line.splitlines():
            if sub:
                yield sub


def iter_rows(
    archive: zipfile.ZipFile,
    entry: zipfile.ZipInfo,
    year: int,
    quarter: int,
    table_kind: str,
    parse_stats: FileParseStats,
) -> Tuple[List[str], Iterator[List[str]]]:
    raw = archive.open(entry)
    text = io.TextIOWrapper(raw, encoding="latin-1", errors="replace", newline="")
    line_iter = iter_fixed_lines(text, year=year, quarter=quarter, table_kind=table_kind)

    try:
        header_line = next(line_iter)
    except StopIteration:
        text.close()
        raw.close()
        return [], iter(())

    header = header_line.rstrip("\r\n").split("$")
    while header and header[-1] == "":
        header.pop()
    header_len = len(header)

    def _row_iter() -> Iterator[List[str]]:
        try:
            for line in line_iter:
                parse_stats.raw_rows += 1
                row = line.rstrip("\r\n").split("$")
                while len(row) > header_len and row[-1] == "":
                    row.pop()
                    parse_stats.trailing_trimmed += 1
                if len(row) < header_len:
                    row.extend([""] * (header_len - len(row)))
                    parse_stats.padded_short += 1
                elif len(row) > header_len:
                    row = row[: header_len - 1] + ["$".join(row[header_len - 1 :])]
                    parse_stats.merged_long += 1
                yield row
        finally:
            text.close()
            raw.close()

    return header, _row_iter()


def find_medascii_dir(meddra_dir: Path) -> Optional[Path]:
    if not meddra_dir.exists():
        return None
    candidate = meddra_dir / "MedAscii"
    if candidate.exists() and candidate.is_dir():
        return candidate
    for d in meddra_dir.rglob("MedAscii"):
        if d.is_dir():
            return d
    return None


def load_meddra_resources(meddra_dir: Optional[Path]) -> MedDRAResources:
    if meddra_dir is None:
        return MedDRAResources(False, {}, {}, {})

    medascii = find_medascii_dir(meddra_dir)
    if medascii is None:
        print(f"[WARN] MedDRA path provided but MedAscii not found: {meddra_dir}", flush=True)
        return MedDRAResources(False, {}, {}, {})

    pt_path = medascii / "pt.asc"
    llt_path = medascii / "llt.asc"
    if not pt_path.exists() or not llt_path.exists():
        print(f"[WARN] MedDRA pt.asc/llt.asc not found in: {medascii}", flush=True)
        return MedDRAResources(False, {}, {}, {})

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
            key = normalize_reaction_term(pt_name)
            if key and key not in pt_name_to_code:
                pt_name_to_code[key] = pt_code
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
            key = normalize_reaction_term(llt_name)
            if key and key not in llt_name_to_pt_code:
                llt_name_to_pt_code[key] = pt_code

    print(
        f"[INFO] MedDRA loaded: PT names={len(pt_name_to_code)}, LLT names={len(llt_name_to_pt_code)}",
        flush=True,
    )
    return MedDRAResources(True, pt_name_to_code, pt_code_to_name, llt_name_to_pt_code)


def standardize_reaction(pt_raw: str, med: MedDRAResources) -> Tuple[str, str, str, str, str]:
    pt_norm = normalize_reaction_term(pt_raw)
    if not pt_norm:
        return "", "", "", "", "EMPTY_PT"

    hint_code = REACTION_MANUAL_CODE_HINTS.get(pt_norm, "")

    if med.loaded:
        pt_code = med.pt_name_to_code.get(pt_norm, "")
        if pt_code:
            pt_name = med.pt_code_to_name.get(pt_code, pt_norm)
            return pt_norm, hint_code, pt_code, pt_name, "PT_EXACT"

        llt_to_pt = med.llt_name_to_pt_code.get(pt_norm, "")
        if llt_to_pt:
            pt_name = med.pt_code_to_name.get(llt_to_pt, "")
            return pt_norm, hint_code, llt_to_pt, pt_name, "LLT_TO_PT"

        if hint_code and hint_code in med.pt_code_to_name:
            pt_name = med.pt_code_to_name[hint_code]
            return pt_norm, hint_code, hint_code, pt_name, "MANUAL_CODE_HINT"

        if hint_code:
            return pt_norm, hint_code, "", "", "MANUAL_HINT_CODE_UNRESOLVED"
        return pt_norm, "", "", "", "UNRESOLVED_WITH_MEDDRA"

    if hint_code:
        return pt_norm, hint_code, "", "", "MANUAL_HINT_NO_MEDDRA"
    return pt_norm, "", "", "", "NO_MEDDRA_DICTIONARY"


def init_db(conn: sqlite3.Connection, rebuild: bool) -> None:
    cur = conn.cursor()
    conn.execute("PRAGMA journal_mode = WAL;")
    conn.execute("PRAGMA synchronous = NORMAL;")
    conn.execute("PRAGMA temp_store = MEMORY;")

    if rebuild:
        cur.executescript(
            """
            DROP TABLE IF EXISTS latest_keys;
            DROP TABLE IF EXISTS case_drug;
            DROP TABLE IF EXISTS case_reaction;
            DROP TABLE IF EXISTS quarter_stats;
            DROP TABLE IF EXISTS reaction_unresolved_top;
            """
        )

    cur.executescript(
        """
        CREATE TABLE IF NOT EXISTS latest_keys (
            year INTEGER NOT NULL,
            quarter INTEGER NOT NULL,
            caseid TEXT NOT NULL,
            primaryid TEXT NOT NULL
        );

        CREATE INDEX IF NOT EXISTS idx_latest_keys_yq ON latest_keys(year, quarter);
        CREATE INDEX IF NOT EXISTS idx_latest_keys_yqpid ON latest_keys(year, quarter, primaryid);

        CREATE TABLE IF NOT EXISTS case_drug (
            year INTEGER NOT NULL,
            quarter INTEGER NOT NULL,
            zip_name TEXT NOT NULL,
            primaryid TEXT NOT NULL,
            caseid TEXT NOT NULL,
            drug_seq TEXT,
            role_cod TEXT,
            drugname_raw TEXT,
            drugname_std TEXT,
            prod_ai_raw TEXT,
            prod_ai_std TEXT,
            drug_key_std TEXT,
            route_raw TEXT,
            route_std TEXT,
            dose_vbm_raw TEXT,
            dose_vbm_std TEXT,
            val_vbm TEXT,
            UNIQUE(primaryid, drug_seq, role_cod, drugname_std, prod_ai_std)
        );

        CREATE INDEX IF NOT EXISTS idx_case_drug_primaryid ON case_drug(primaryid);
        CREATE INDEX IF NOT EXISTS idx_case_drug_drugkey ON case_drug(drug_key_std);

        CREATE TABLE IF NOT EXISTS case_reaction (
            year INTEGER NOT NULL,
            quarter INTEGER NOT NULL,
            zip_name TEXT NOT NULL,
            primaryid TEXT NOT NULL,
            caseid TEXT NOT NULL,
            pt_raw TEXT,
            pt_norm TEXT,
            pt_manual_code_hint TEXT,
            pt_meddra_code TEXT,
            pt_meddra_pt TEXT,
            pt_mapping_source TEXT,
            drug_rec_act_raw TEXT,
            drug_rec_act_std TEXT,
            UNIQUE(primaryid, pt_norm)
        );

        CREATE INDEX IF NOT EXISTS idx_case_reaction_primaryid ON case_reaction(primaryid);
        CREATE INDEX IF NOT EXISTS idx_case_reaction_ptnorm ON case_reaction(pt_norm);

        CREATE TABLE IF NOT EXISTS quarter_stats (
            table_kind TEXT NOT NULL,
            zip_name TEXT NOT NULL,
            year INTEGER NOT NULL,
            quarter INTEGER NOT NULL,
            entry_name TEXT,
            raw_rows INTEGER,
            filtered_rows INTEGER,
            inserted_rows INTEGER,
            skipped_no_primaryid INTEGER,
            trailing_trimmed INTEGER,
            padded_short INTEGER,
            merged_long INTEGER,
            caseid_overridden INTEGER
        );

        CREATE TABLE IF NOT EXISTS reaction_unresolved_top (
            pt_norm TEXT NOT NULL,
            cnt INTEGER NOT NULL
        );
        """
    )
    conn.commit()


def load_latest_keys(conn: sqlite3.Connection, key_csv: Path, rebuild: bool) -> None:
    cur = conn.cursor()
    if rebuild:
        cur.execute("DELETE FROM latest_keys")
        conn.commit()

    existing = cur.execute("SELECT COUNT(*) FROM latest_keys").fetchone()[0]
    if existing > 0:
        print(f"[INFO] latest_keys already loaded: {existing}", flush=True)
        return

    with key_csv.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        batch: List[Tuple[int, int, str, str]] = []
        for row in reader:
            try:
                year = int(row["year"])
                quarter = int(row["quarter"])
            except Exception:
                continue
            caseid = row.get("caseid", "").strip()
            primaryid = row.get("primaryid", "").strip()
            if not caseid or not primaryid:
                continue
            batch.append((year, quarter, caseid, primaryid))
            if len(batch) >= 200000:
                conn.executemany(
                    "INSERT INTO latest_keys(year,quarter,caseid,primaryid) VALUES (?,?,?,?)",
                    batch,
                )
                conn.commit()
                batch.clear()
        if batch:
            conn.executemany(
                "INSERT INTO latest_keys(year,quarter,caseid,primaryid) VALUES (?,?,?,?)",
                batch,
            )
            conn.commit()

    total = cur.execute("SELECT COUNT(*) FROM latest_keys").fetchone()[0]
    print(f"[INFO] latest_keys loaded: {total}", flush=True)


def fetch_selected_map(conn: sqlite3.Connection, year: int, quarter: int) -> Dict[str, str]:
    cur = conn.cursor()
    rows = cur.execute(
        "SELECT primaryid, caseid FROM latest_keys WHERE year=? AND quarter=?",
        (year, quarter),
    )
    return {r[0]: r[1] for r in rows}


def ingest_drug_file(
    conn: sqlite3.Connection,
    archive: zipfile.ZipFile,
    entry: zipfile.ZipInfo,
    info: ZipInfo,
    selected_map: Dict[str, str],
    batch_size: int = 50000,
) -> FileParseStats:
    st = FileParseStats()
    header, rows = iter_rows(
        archive=archive,
        entry=entry,
        year=info.year,
        quarter=info.quarter,
        table_kind="drug",
        parse_stats=st,
    )
    if not header:
        return st

    idx = build_header_index(header, DRUG_SYNONYMS)
    insert_sql = (
        "INSERT OR IGNORE INTO case_drug("
        "year,quarter,zip_name,primaryid,caseid,drug_seq,role_cod,drugname_raw,drugname_std,"
        "prod_ai_raw,prod_ai_std,drug_key_std,route_raw,route_std,dose_vbm_raw,dose_vbm_std,val_vbm"
        ") VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
    )

    batch: List[Tuple[object, ...]] = []
    before_changes = conn.total_changes

    for row in rows:
        primaryid = cell(row, idx["primaryid"]).strip()
        if not primaryid:
            st.skipped_no_primaryid += 1
            continue
        selected_caseid = selected_map.get(primaryid)
        if selected_caseid is None:
            continue

        row_caseid = cell(row, idx["caseid"]).strip()
        caseid = selected_caseid
        if row_caseid and row_caseid != selected_caseid:
            st.caseid_overridden += 1

        drug_seq = normalize_text(cell(row, idx["drug_seq"]), upper=False)
        role_cod = normalize_text(cell(row, idx["role_cod"]), upper=True)

        drugname_raw = cell(row, idx["drugname"])
        prod_ai_raw = cell(row, idx["prod_ai"])
        route_raw = cell(row, idx["route"])
        dose_vbm_raw = cell(row, idx["dose_vbm"])
        val_vbm = normalize_text(cell(row, idx["val_vbm"]), upper=False)

        drugname_std = normalize_drug_term(drugname_raw)
        prod_ai_std = normalize_drug_term(prod_ai_raw)
        route_std = normalize_text(route_raw, upper=True)
        dose_vbm_std = normalize_text(dose_vbm_raw, upper=True)

        drug_key_std = prod_ai_std if prod_ai_std else drugname_std

        st.filtered_rows += 1
        batch.append(
            (
                info.year,
                info.quarter,
                info.path.name,
                primaryid,
                caseid,
                drug_seq,
                role_cod,
                drugname_raw,
                drugname_std,
                prod_ai_raw,
                prod_ai_std,
                drug_key_std,
                route_raw,
                route_std,
                dose_vbm_raw,
                dose_vbm_std,
                val_vbm,
            )
        )

        if len(batch) >= batch_size:
            conn.executemany(insert_sql, batch)
            conn.commit()
            batch.clear()

    if batch:
        conn.executemany(insert_sql, batch)
        conn.commit()

    st.inserted_rows = conn.total_changes - before_changes
    return st


def ingest_reac_file(
    conn: sqlite3.Connection,
    archive: zipfile.ZipFile,
    entry: zipfile.ZipInfo,
    info: ZipInfo,
    selected_map: Dict[str, str],
    med: MedDRAResources,
    unresolved_counter: Counter,
    mapping_source_counter: Counter,
    batch_size: int = 50000,
) -> FileParseStats:
    st = FileParseStats()
    header, rows = iter_rows(
        archive=archive,
        entry=entry,
        year=info.year,
        quarter=info.quarter,
        table_kind="reac",
        parse_stats=st,
    )
    if not header:
        return st

    idx = build_header_index(header, REAC_SYNONYMS)
    insert_sql = (
        "INSERT OR IGNORE INTO case_reaction("
        "year,quarter,zip_name,primaryid,caseid,pt_raw,pt_norm,pt_manual_code_hint,"
        "pt_meddra_code,pt_meddra_pt,pt_mapping_source,drug_rec_act_raw,drug_rec_act_std"
        ") VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)"
    )

    batch: List[Tuple[object, ...]] = []
    before_changes = conn.total_changes

    for row in rows:
        primaryid = cell(row, idx["primaryid"]).strip()
        if not primaryid:
            st.skipped_no_primaryid += 1
            continue

        selected_caseid = selected_map.get(primaryid)
        if selected_caseid is None:
            continue

        row_caseid = cell(row, idx["caseid"]).strip()
        caseid = selected_caseid
        if row_caseid and row_caseid != selected_caseid:
            st.caseid_overridden += 1

        pt_raw = cell(row, idx["pt"])
        pt_norm, code_hint, meddra_code, meddra_pt, map_source = standardize_reaction(pt_raw, med)
        mapping_source_counter[map_source] += 1
        if map_source in {"NO_MEDDRA_DICTIONARY", "UNRESOLVED_WITH_MEDDRA", "MANUAL_HINT_CODE_UNRESOLVED"}:
            if pt_norm:
                unresolved_counter[pt_norm] += 1

        drug_rec_act_raw = cell(row, idx["drug_rec_act"])
        drug_rec_act_std = normalize_text(drug_rec_act_raw, upper=True)

        st.filtered_rows += 1
        batch.append(
            (
                info.year,
                info.quarter,
                info.path.name,
                primaryid,
                caseid,
                pt_raw,
                pt_norm,
                code_hint,
                meddra_code,
                meddra_pt,
                map_source,
                drug_rec_act_raw,
                drug_rec_act_std,
            )
        )

        if len(batch) >= batch_size:
            conn.executemany(insert_sql, batch)
            conn.commit()
            batch.clear()

    if batch:
        conn.executemany(insert_sql, batch)
        conn.commit()

    st.inserted_rows = conn.total_changes - before_changes
    return st


def insert_quarter_stats(
    conn: sqlite3.Connection,
    table_kind: str,
    info: ZipInfo,
    entry_name: str,
    st: FileParseStats,
) -> None:
    conn.execute(
        """
        INSERT INTO quarter_stats(
            table_kind, zip_name, year, quarter, entry_name,
            raw_rows, filtered_rows, inserted_rows, skipped_no_primaryid,
            trailing_trimmed, padded_short, merged_long, caseid_overridden
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)
        """,
        (
            table_kind,
            info.path.name,
            info.year,
            info.quarter,
            entry_name,
            st.raw_rows,
            st.filtered_rows,
            st.inserted_rows,
            st.skipped_no_primaryid,
            st.trailing_trimmed,
            st.padded_short,
            st.merged_long,
            st.caseid_overridden,
        ),
    )
    conn.commit()


def already_processed(conn: sqlite3.Connection, table_kind: str, info: ZipInfo) -> bool:
    row = conn.execute(
        """
        SELECT COUNT(*)
        FROM quarter_stats
        WHERE table_kind = ? AND zip_name = ? AND year = ? AND quarter = ?
        """,
        (table_kind, info.path.name, info.year, info.quarter),
    ).fetchone()
    return bool(row and row[0] > 0)


def export_query_to_csv(conn: sqlite3.Connection, sql: str, out_path: Path) -> int:
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
    p = argparse.ArgumentParser(description="Step2 build standardized case-drug/reaction")
    p.add_argument(
        "--data-root",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "raw_data" / "faers_quarterly_archives",
    )
    p.add_argument(
        "--step1-keys",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "outputs" / "step1" / "faers_latest_case_keys.csv",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "outputs" / "step2",
    )
    p.add_argument(
        "--meddra-dir",
        type=Path,
        default=None,
        help="Optional path to MedDRA directory containing MedAscii/pt.asc and llt.asc",
    )
    p.add_argument("--rebuild", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if not args.data_root.exists():
        print(f"[ERROR] data-root not found: {args.data_root}", file=sys.stderr)
        return 1
    if not args.step1_keys.exists():
        print(f"[ERROR] step1 keys not found: {args.step1_keys}", file=sys.stderr)
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    db_path = args.output_dir / "faers_step2.sqlite"

    started = dt.datetime.now().isoformat(timespec="seconds")

    med = load_meddra_resources(args.meddra_dir)

    conn = sqlite3.connect(db_path)
    unresolved_counter: Counter = Counter()
    mapping_source_counter: Counter = Counter()

    try:
        init_db(conn, rebuild=args.rebuild)
        load_latest_keys(conn, args.step1_keys, rebuild=args.rebuild)

        zip_infos = collect_zip_files(args.data_root)
        if not zip_infos:
            print("[ERROR] No zip files found", file=sys.stderr)
            return 1

        print(f"[INFO] Processing {len(zip_infos)} zip files...", flush=True)

        for info in zip_infos:
            selected_map = fetch_selected_map(conn, info.year, info.quarter)
            if not selected_map:
                continue

            with zipfile.ZipFile(info.path, "r") as arc:
                drug_entry = find_entry(arc, DRUG_ENTRY_RE)
                reac_entry = find_entry(arc, REAC_ENTRY_RE)

                if already_processed(conn, "drug", info):
                    pass
                elif drug_entry is not None:
                    drug_stats = ingest_drug_file(
                        conn=conn,
                        archive=arc,
                        entry=drug_entry,
                        info=info,
                        selected_map=selected_map,
                    )
                    insert_quarter_stats(conn, "drug", info, drug_entry.filename, drug_stats)
                else:
                    insert_quarter_stats(conn, "drug", info, "", FileParseStats())

                if already_processed(conn, "reac", info):
                    pass
                elif reac_entry is not None:
                    reac_stats = ingest_reac_file(
                        conn=conn,
                        archive=arc,
                        entry=reac_entry,
                        info=info,
                        selected_map=selected_map,
                        med=med,
                        unresolved_counter=unresolved_counter,
                        mapping_source_counter=mapping_source_counter,
                    )
                    insert_quarter_stats(conn, "reac", info, reac_entry.filename, reac_stats)
                else:
                    insert_quarter_stats(conn, "reac", info, "", FileParseStats())

            print(
                f"[Q] {info.path.name}: keys={len(selected_map)} processed",
                flush=True,
            )

        # Persist top unresolved reactions
        conn.execute("DELETE FROM reaction_unresolved_top")
        top_unresolved = unresolved_counter.most_common(2000)
        if top_unresolved:
            conn.executemany(
                "INSERT INTO reaction_unresolved_top(pt_norm,cnt) VALUES (?,?)",
                top_unresolved,
            )
            conn.commit()

        # Export outputs
        case_drug_csv = args.output_dir / "case_drug_filtered_std.csv"
        case_reac_csv = args.output_dir / "case_reaction_filtered_std.csv"
        quarter_stats_csv = args.output_dir / "step2_quarter_stats.csv"
        unresolved_csv = args.output_dir / "reaction_unresolved_top.csv"

        n_drug = export_query_to_csv(
            conn,
            "SELECT * FROM case_drug ORDER BY year, quarter, primaryid",
            case_drug_csv,
        )
        n_reac = export_query_to_csv(
            conn,
            "SELECT * FROM case_reaction ORDER BY year, quarter, primaryid",
            case_reac_csv,
        )
        export_query_to_csv(
            conn,
            "SELECT * FROM quarter_stats ORDER BY year, quarter, table_kind",
            quarter_stats_csv,
        )
        export_query_to_csv(
            conn,
            "SELECT * FROM reaction_unresolved_top ORDER BY cnt DESC, pt_norm",
            unresolved_csv,
        )

        finished = dt.datetime.now().isoformat(timespec="seconds")
        summary = {
            "started_at": started,
            "finished_at": finished,
            "meddra_loaded": med.loaded,
            "zip_files_total": len(zip_infos),
            "latest_keys_rows": scalar(conn, "SELECT COUNT(*) FROM latest_keys"),
            "case_drug_rows": n_drug,
            "case_reaction_rows": n_reac,
            "case_drug_distinct_primaryid": scalar(conn, "SELECT COUNT(DISTINCT primaryid) FROM case_drug"),
            "case_reaction_distinct_primaryid": scalar(conn, "SELECT COUNT(DISTINCT primaryid) FROM case_reaction"),
            "reaction_mapping_source_counts": dict(mapping_source_counter),
            "reaction_unresolved_unique_terms": len(unresolved_counter),
            "reaction_unresolved_total_rows": int(sum(unresolved_counter.values())),
        }

        (args.output_dir / "step2_report.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

        print("[DONE] Step2 completed", flush=True)
        print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    raise SystemExit(main())

