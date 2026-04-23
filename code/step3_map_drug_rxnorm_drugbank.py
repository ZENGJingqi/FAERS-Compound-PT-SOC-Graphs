#!/usr/bin/env python3
"""Step 3: Normalize drug terms and link them to RxNorm and DrugBank.

Purpose:
- collapse standardized Step2 drug terms to unique normalized drug strings
- map drug names to RxNorm-compatible ingredients or products
- connect mapped drug names to DrugBank identifiers and InChI records

Inputs:
- outputs/step2/faers_step2.sqlite
- reference_data/RxNorm_full_03022026/
- reference_data/drugbank_5.1.15/

Outputs:
- outputs/step3/faers_step3.sqlite
- outputs/step3/*.csv
- outputs/step3/*.json
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import re
import sqlite3
import unicodedata
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

SAB_PRIORITY = {
    "RXNORM": 1,
    "MTHSPL": 2,
    "MMSL": 3,
    "MMX": 4,
    "GS": 5,
    "VANDF": 6,
    "NDDF": 7,
    "ATC": 8,
    "DRUGBANK": 9,
    "USP": 10,
    "CVX": 11,
}

TTY_PRIORITY = {
    "IN": 1,
    "PIN": 2,
    "MIN": 3,
    "SCD": 4,
    "SBD": 5,
    "SCDC": 6,
    "SBDC": 7,
    "GPCK": 8,
    "BPCK": 9,
    "SCDG": 10,
    "SBDG": 11,
    "BN": 12,
    "DF": 13,
}

STRENGTH_RE = re.compile(
    r"\b\d+(?:\.\d+)?\s*(?:MG|MCG|UG|G|GRAM|GRAMS|KG|ML|MILLILITER|MILLILITERS|MILLILITRE|MILLILITRES|IU|UNITS|%)"
    r"(?:\s*/\s*\d+(?:\.\d+)?\s*(?:MG|MCG|UG|G|GRAM|GRAMS|KG|ML|MILLILITER|MILLILITERS|MILLILITRE|MILLILITRES|IU|UNITS|%))?\b"
)
FORM_WORD_RE = re.compile(
    r"\b(?:TABLET|TABLETS|CAPSULE|CAPSULES|INJECTION|INJECTABLE|SOLUTION|SUSPENSION|SYRUP|CREAM|OINTMENT|"
    r"UNKNOWN|UNK|NOS|FORMULATION|GENERIC|BLINDED|DOSE|DOSES|MG|MCG|ML|IU)\b"
)


@dataclass
class TermStat:
    term_std: str
    term_norm: str
    n_rows: int
    n_cases: int


@dataclass(frozen=True)
class RxCandidate:
    rxcui: str
    sab: str
    tty: str
    code: str
    string: str


@dataclass(frozen=True)
class ManualCorrection:
    action: str
    target_term: str
    note: str


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
    t = t.replace("\\", " / ")
    t = t.strip("\"'")
    t = re.sub(r"\s*([,;/+&()\[\]-])\s*", r" \1 ", t)
    t = re.sub(r"\s+", " ", t).strip()
    return t


def clean_drug_term(s: str) -> str:
    t = normalize_drug_term(s)
    if not t:
        return ""
    t = STRENGTH_RE.sub(" ", t)
    t = FORM_WORD_RE.sub(" ", t)
    t = re.sub(r"\b(?:\d+(?:\.\d+)?)\b", " ", t)
    t = re.sub(r"\s*([,;/+()\[\]-])\s*", " ", t)
    t = re.sub(r"\s+", " ", t).strip()
    return t


def build_variants(term_std: str) -> List[Tuple[int, str, str]]:
    variants: List[Tuple[int, str, str]] = []
    seen: Set[str] = set()

    def add(pri: int, stage: str, value: str) -> None:
        v = normalize_drug_term(value)
        if not v:
            return
        if v in seen:
            return
        seen.add(v)
        variants.append((pri, stage, v))

    base = normalize_drug_term(term_std)
    add(1, "EXACT", base)

    cleaned = clean_drug_term(term_std)
    add(2, "CLEANED", cleaned)

    for m in re.finditer(r"\(([^()]{2,})\)", base):
        add(3, "PAREN", m.group(1))

    for sep in ["/", "+", ",", ";", "&"]:
        if sep in base:
            for part in base.split(sep):
                add(4, "SPLIT", part)

    return sorted(variants, key=lambda x: (x[0], x[1], x[2]))


def candidate_rank(c: RxCandidate, target_norm: str) -> Tuple[int, int, int, int, int]:
    sab_p = SAB_PRIORITY.get(c.sab, 99)
    tty_p = TTY_PRIORITY.get(c.tty, 99)
    s_norm = normalize_drug_term(c.string)
    exact = 0 if s_norm == target_norm else 1
    len_diff = abs(len(s_norm) - len(target_norm))
    return (sab_p, tty_p, exact, len_diff, len(s_norm))


def init_db(conn: sqlite3.Connection, rebuild: bool) -> None:
    cur = conn.cursor()
    conn.execute("PRAGMA journal_mode = WAL;")
    conn.execute("PRAGMA synchronous = NORMAL;")
    conn.execute("PRAGMA temp_store = FILE;")

    if rebuild:
        cur.executescript(
            """
            DROP TABLE IF EXISTS drug_term_stats;
            DROP TABLE IF EXISTS drug_term_final;
            DROP TABLE IF EXISTS unmapped_top;
            DROP TABLE IF EXISTS term_stage_counts;
            """
        )

    cur.executescript(
        """
        CREATE TABLE IF NOT EXISTS drug_term_stats (
            term_std TEXT PRIMARY KEY,
            term_norm TEXT NOT NULL,
            n_rows INTEGER NOT NULL,
            n_cases INTEGER NOT NULL
        );

        CREATE TABLE IF NOT EXISTS drug_term_final (
            term_std TEXT PRIMARY KEY,
            term_norm TEXT NOT NULL,
            n_rows INTEGER NOT NULL,
            n_cases INTEGER NOT NULL,
            manual_action TEXT,
            manual_target TEXT,
            mapping_stage TEXT,
            candidate_count INTEGER,
            rxcui TEXT,
            sab TEXT,
            tty TEXT,
            code TEXT,
            rx_string TEXT,
            dbid TEXT,
            drugbank_name TEXT,
            inchi TEXT
        );

        CREATE TABLE IF NOT EXISTS unmapped_top (
            term_std TEXT,
            term_norm TEXT,
            n_rows INTEGER,
            n_cases INTEGER
        );

        CREATE TABLE IF NOT EXISTS term_stage_counts (
            mapping_stage TEXT PRIMARY KEY,
            term_count INTEGER,
            row_count INTEGER
        );
        """
    )
    conn.commit()


def load_term_stats_from_step2(step2_db: Path) -> List[TermStat]:
    conn = sqlite3.connect(f"file:{step2_db}?mode=ro", uri=True)
    cur = conn.cursor()
    q = (
        "SELECT drug_key_std, COUNT(*) AS n_rows, COUNT(DISTINCT primaryid) AS n_cases "
        "FROM case_drug WHERE TRIM(drug_key_std) <> '' GROUP BY drug_key_std"
    )
    out: List[TermStat] = []
    for term_std, n_rows, n_cases in cur.execute(q):
        tnorm = normalize_drug_term(term_std)
        if not tnorm:
            continue
        out.append(TermStat(term_std=term_std, term_norm=tnorm, n_rows=int(n_rows), n_cases=int(n_cases)))
    conn.close()
    return out


def save_term_stats(conn: sqlite3.Connection, terms: Sequence[TermStat]) -> None:
    conn.executemany(
        "INSERT OR REPLACE INTO drug_term_stats(term_std,term_norm,n_rows,n_cases) VALUES (?,?,?,?)",
        [(t.term_std, t.term_norm, t.n_rows, t.n_cases) for t in terms],
    )
    conn.commit()


def load_drugbank_maps(drugbank_dir: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    dname_path = drugbank_dir / "XML_dbid_dname.csv"
    inchi_path = drugbank_dir / "XML_dbid_inchi.csv"

    dbname: Dict[str, str] = {}
    dbinchi: Dict[str, str] = {}

    with dname_path.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            dbid = normalize_text(row.get("DrugID", ""), upper=True)
            name = normalize_text(row.get("DrugName", ""), upper=False)
            if dbid and dbid not in dbname:
                dbname[dbid] = name

    with inchi_path.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            dbid = normalize_text(row.get("DrugID", ""), upper=True)
            inchi = normalize_text(row.get("InChI", ""), upper=False)
            if dbid and inchi and dbid not in dbinchi:
                dbinchi[dbid] = inchi

    return dbname, dbinchi


def build_drugbank_name_index(dbname: Dict[str, str]) -> Dict[str, Set[str]]:
    out: Dict[str, Set[str]] = defaultdict(set)
    for dbid, name in dbname.items():
        norm = normalize_drug_term(name)
        if norm:
            out[norm].add(dbid)
        cleaned = clean_drug_term(name)
        if cleaned:
            out[cleaned].add(dbid)
    return out


def load_manual_corrections(path: Optional[Path]) -> Dict[str, ManualCorrection]:
    if path is None or not path.exists():
        return {}

    out: Dict[str, ManualCorrection] = {}
    with path.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            term_raw = row.get("term_std", "")
            action_raw = row.get("action", "")
            target_raw = row.get("target_term", "")
            note = (row.get("note", "") or "").strip()

            key = normalize_drug_term(term_raw)
            if not key:
                continue

            action = normalize_text(action_raw, upper=True)
            if not action:
                continue
            if action not in {"MAP_TO_TERM", "EXCLUDE_NON_DRUG"}:
                continue

            target = normalize_drug_term(target_raw) if target_raw else ""
            out[key] = ManualCorrection(action=action, target_term=target, note=note)
    return out


def choose_dbid_by_name_match(
    candidate_terms: Iterable[str],
    dbname_index: Dict[str, Set[str]],
    dbname: Dict[str, str],
    dbinchi: Dict[str, str],
) -> Tuple[str, str, str]:
    normalized_terms: List[str] = []
    seen: Set[str] = set()
    for t in candidate_terms:
        n = normalize_drug_term(t)
        c = clean_drug_term(t)
        for x in (n, c):
            if x and x not in seen:
                seen.add(x)
                normalized_terms.append(x)

    for term in normalized_terms:
        dbids = sorted(dbname_index.get(term, set()))
        for dbid in dbids:
            inchi = dbinchi.get(dbid, "")
            if inchi:
                return dbid, dbname.get(dbid, ""), inchi
    return "", "", ""


def build_rxnorm_candidate_maps(
    rxnconso_path: Path,
    variant_keys: Set[str],
) -> Tuple[Dict[str, List[RxCandidate]], Dict[str, Set[str]], Dict[str, int]]:
    rx_by_term: Dict[str, List[RxCandidate]] = defaultdict(list)
    rxcui_to_dbids: Dict[str, Set[str]] = defaultdict(set)

    scan_stats = {
        "lines_total": 0,
        "lines_kept_variant": 0,
        "lines_english_active": 0,
        "drugbank_links": 0,
    }

    with rxnconso_path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        for line in f:
            scan_stats["lines_total"] += 1
            parts = line.rstrip("\n\r").split("|")
            if len(parts) < 17:
                continue
            rxcui = parts[0].strip()
            lat = parts[1].strip()
            sab = parts[11].strip()
            tty = parts[12].strip()
            code = parts[13].strip()
            string = parts[14].strip()
            suppress = parts[16].strip()

            if sab == "DRUGBANK" and code.startswith("DB") and rxcui:
                rxcui_to_dbids[rxcui].add(code)
                scan_stats["drugbank_links"] += 1

            if lat != "ENG" or suppress != "N":
                continue
            scan_stats["lines_english_active"] += 1

            term_norm = normalize_drug_term(string)
            if not term_norm or term_norm not in variant_keys:
                continue

            cand = RxCandidate(rxcui=rxcui, sab=sab, tty=tty, code=code, string=string)
            rx_by_term[term_norm].append(cand)
            scan_stats["lines_kept_variant"] += 1

    # de-duplicate candidates per term
    for key, vals in list(rx_by_term.items()):
        seen = set()
        uniq: List[RxCandidate] = []
        for c in vals:
            k = (c.rxcui, c.sab, c.tty, c.code, c.string)
            if k in seen:
                continue
            seen.add(k)
            uniq.append(c)
        rx_by_term[key] = uniq

    return rx_by_term, rxcui_to_dbids, scan_stats


def choose_dbid_for_rxcui(
    rxcui: str,
    rxcui_to_dbids: Dict[str, Set[str]],
    dbname: Dict[str, str],
    dbinchi: Dict[str, str],
) -> Tuple[str, str, str]:
    dbids = sorted(rxcui_to_dbids.get(rxcui, set()))
    if not dbids:
        return "", "", ""

    # Prefer entries with InChI
    with_inchi = [d for d in dbids if d in dbinchi]
    chosen = with_inchi[0] if with_inchi else dbids[0]
    return chosen, dbname.get(chosen, ""), dbinchi.get(chosen, "")


def map_terms(
    terms: Sequence[TermStat],
    rx_by_term: Dict[str, List[RxCandidate]],
    rxcui_to_dbids: Dict[str, Set[str]],
    dbname: Dict[str, str],
    dbinchi: Dict[str, str],
    dbname_index: Dict[str, Set[str]],
    manual_corrections: Dict[str, ManualCorrection],
) -> Tuple[List[Tuple], Dict[str, int], Dict[str, int]]:
    final_rows: List[Tuple] = []
    stage_term_count: Counter = Counter()
    stage_row_count: Counter = Counter()

    for t in terms:
        term_key = normalize_drug_term(t.term_std)
        corr = manual_corrections.get(term_key)
        manual_action = ""
        manual_target = ""
        base_term = t.term_std

        if corr:
            manual_action = corr.action
            manual_target = corr.target_term
            if corr.action == "EXCLUDE_NON_DRUG":
                stage = "MANUAL_EXCLUDE"
                final_rows.append(
                    (
                        t.term_std,
                        t.term_norm,
                        t.n_rows,
                        t.n_cases,
                        manual_action,
                        manual_target,
                        stage,
                        0,
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                        "",
                    )
                )
                stage_term_count[stage] += 1
                stage_row_count[stage] += t.n_rows
                continue
            if corr.action == "MAP_TO_TERM" and corr.target_term:
                base_term = corr.target_term

        variants = build_variants(base_term)

        chosen_stage = "UNMAPPED"
        stage_candidates: List[RxCandidate] = []

        grouped: Dict[int, List[Tuple[str, str]]] = defaultdict(list)
        for pri, stage, norm_v in variants:
            grouped[pri].append((stage, norm_v))

        for pri in sorted(grouped.keys()):
            pool: List[RxCandidate] = []
            stage_name = grouped[pri][0][0]
            for _, v in grouped[pri]:
                pool.extend(rx_by_term.get(v, []))
            if pool:
                # unique in stage
                uniq = []
                seen = set()
                for c in pool:
                    k = (c.rxcui, c.sab, c.tty, c.code, c.string)
                    if k in seen:
                        continue
                    seen.add(k)
                    uniq.append(c)
                stage_candidates = uniq
                chosen_stage = stage_name
                break

        if not stage_candidates:
            dbid, dbnm, inchi = choose_dbid_by_name_match(
                candidate_terms=[base_term, t.term_std, t.term_norm],
                dbname_index=dbname_index,
                dbname=dbname,
                dbinchi=dbinchi,
            )
            if dbid:
                stage = "DB_NAME_EXACT"
                if manual_action == "MAP_TO_TERM":
                    stage = f"MANUAL_{stage}"
                final_rows.append(
                    (
                        t.term_std,
                        t.term_norm,
                        t.n_rows,
                        t.n_cases,
                        manual_action,
                        manual_target,
                        stage,
                        1,
                        "",
                        "",
                        "",
                        "",
                        "",
                        dbid,
                        dbnm,
                        inchi,
                    )
                )
                stage_term_count[stage] += 1
                stage_row_count[stage] += t.n_rows
                continue

            stage = "UNMAPPED"
            if manual_action == "MAP_TO_TERM":
                stage = "MANUAL_UNMAPPED"
            final_rows.append(
                (
                    t.term_std,
                    t.term_norm,
                    t.n_rows,
                    t.n_cases,
                    manual_action,
                    manual_target,
                    stage,
                    0,
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                )
            )
            stage_term_count[stage] += 1
            stage_row_count[stage] += t.n_rows
            continue

        best = sorted(
            stage_candidates,
            key=lambda c: candidate_rank(c, t.term_norm),
        )[0]

        dbid, dbnm, inchi = choose_dbid_for_rxcui(best.rxcui, rxcui_to_dbids, dbname, dbinchi)
        stage = chosen_stage
        if manual_action == "MAP_TO_TERM":
            stage = f"MANUAL_{stage}"

        final_rows.append(
            (
                t.term_std,
                t.term_norm,
                t.n_rows,
                t.n_cases,
                manual_action,
                manual_target,
                stage,
                len(stage_candidates),
                best.rxcui,
                best.sab,
                best.tty,
                best.code,
                best.string,
                dbid,
                dbnm,
                inchi,
            )
        )
        stage_term_count[stage] += 1
        stage_row_count[stage] += t.n_rows

    return final_rows, dict(stage_term_count), dict(stage_row_count)


def save_results(
    conn: sqlite3.Connection,
    final_rows: Sequence[Tuple],
    stage_term_count: Dict[str, int],
    stage_row_count: Dict[str, int],
) -> None:
    conn.execute("DELETE FROM drug_term_final")
    conn.execute("DELETE FROM unmapped_top")
    conn.execute("DELETE FROM term_stage_counts")

    conn.executemany(
        """
        INSERT INTO drug_term_final(
            term_std,term_norm,n_rows,n_cases,manual_action,manual_target,mapping_stage,candidate_count,
            rxcui,sab,tty,code,rx_string,dbid,drugbank_name,inchi
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """,
        final_rows,
    )

    conn.executemany(
        "INSERT INTO term_stage_counts(mapping_stage,term_count,row_count) VALUES (?,?,?)",
        [(k, int(stage_term_count.get(k, 0)), int(stage_row_count.get(k, 0))) for k in sorted(stage_term_count.keys())],
    )

    conn.execute(
        """
        INSERT INTO unmapped_top(term_std,term_norm,n_rows,n_cases)
        SELECT term_std,term_norm,n_rows,n_cases
        FROM drug_term_final
        WHERE COALESCE(inchi,'') = ''
        ORDER BY n_rows DESC
        LIMIT 5000
        """
    )
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
    p = argparse.ArgumentParser(description="Step3 map drug terms to RxNorm and DrugBank")
    p.add_argument(
        "--step2-db",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "outputs" / "step2" / "faers_step2.sqlite",
    )
    p.add_argument(
        "--rxnconso",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "reference_data" / "RxNorm_full_03022026" / "rrf" / "RXNCONSO.RRF",
    )
    p.add_argument(
        "--drugbank-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "reference_data" / "drugbank_5.1.15",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "outputs" / "step3",
    )
    p.add_argument(
        "--manual-corrections",
        type=Path,
        default=Path(
            str(Path(__file__).resolve().parents[1] / "resources" / "step3_manual" / "corrections_batch1_top300.csv")
        ),
    )
    p.add_argument("--rebuild", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if not args.step2_db.exists():
        print(f"[ERROR] step2 db not found: {args.step2_db}")
        return 1
    if not args.rxnconso.exists():
        print(f"[ERROR] RXNCONSO not found: {args.rxnconso}")
        return 1
    if not args.drugbank_dir.exists():
        print(f"[ERROR] drugbank dir not found: {args.drugbank_dir}")
        return 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    step3_db = args.output_dir / "faers_step3.sqlite"
    started = dt.datetime.now().isoformat(timespec="seconds")

    conn = sqlite3.connect(step3_db)
    try:
        init_db(conn, rebuild=args.rebuild)

        print("[INFO] loading unique drug terms from step2 ...", flush=True)
        terms = load_term_stats_from_step2(args.step2_db)
        save_term_stats(conn, terms)
        print(f"[INFO] terms loaded: {len(terms)}", flush=True)

        manual_corrections = load_manual_corrections(args.manual_corrections)
        print(f"[INFO] manual corrections loaded: {len(manual_corrections)}", flush=True)

        print("[INFO] building term variants ...", flush=True)
        variant_keys: Set[str] = set()
        for t in terms:
            term_key = normalize_drug_term(t.term_std)
            corr = manual_corrections.get(term_key)
            base_term = corr.target_term if (corr and corr.action == "MAP_TO_TERM" and corr.target_term) else t.term_std
            for _, _, v in build_variants(base_term):
                variant_keys.add(v)
        print(f"[INFO] unique variant keys: {len(variant_keys)}", flush=True)

        print("[INFO] loading DrugBank maps ...", flush=True)
        dbname, dbinchi = load_drugbank_maps(args.drugbank_dir)
        dbname_index = build_drugbank_name_index(dbname)
        print(
            f"[INFO] DrugBank names={len(dbname)} inchi={len(dbinchi)} name_keys={len(dbname_index)}",
            flush=True,
        )

        print("[INFO] scanning RXNCONSO and building candidate maps ...", flush=True)
        rx_by_term, rxcui_to_dbids, scan_stats = build_rxnorm_candidate_maps(
            args.rxnconso,
            variant_keys,
        )
        print(
            f"[INFO] RX candidates terms={len(rx_by_term)} rxcui->dbid={len(rxcui_to_dbids)}",
            flush=True,
        )

        print("[INFO] mapping terms ...", flush=True)
        final_rows, stage_term_count, stage_row_count = map_terms(
            terms,
            rx_by_term,
            rxcui_to_dbids,
            dbname,
            dbinchi,
            dbname_index,
            manual_corrections,
        )

        print("[INFO] saving results ...", flush=True)
        save_results(conn, final_rows, stage_term_count, stage_row_count)

        # Export
        n_final = export_csv(
            conn,
            """
            SELECT term_std,term_norm,n_rows,n_cases,manual_action,manual_target,mapping_stage,candidate_count,
                   rxcui,sab,tty,code,rx_string,dbid,drugbank_name,inchi
            FROM drug_term_final
            ORDER BY n_rows DESC, term_std
            """,
            args.output_dir / "drug_term_final.csv",
        )
        export_csv(
            conn,
            "SELECT * FROM unmapped_top ORDER BY n_rows DESC, term_std",
            args.output_dir / "unmapped_drug_terms_top.csv",
        )
        export_csv(
            conn,
            "SELECT * FROM term_stage_counts ORDER BY row_count DESC, mapping_stage",
            args.output_dir / "term_stage_counts.csv",
        )

        # Coverage summary
        term_total = scalar(conn, "SELECT COUNT(*) FROM drug_term_final")
        term_mapped_rxcui = scalar(conn, "SELECT COUNT(*) FROM drug_term_final WHERE COALESCE(rxcui,'') <> ''")
        term_mapped_dbid = scalar(conn, "SELECT COUNT(*) FROM drug_term_final WHERE COALESCE(dbid,'') <> ''")
        term_mapped_inchi = scalar(conn, "SELECT COUNT(*) FROM drug_term_final WHERE COALESCE(inchi,'') <> ''")
        term_manual_excluded = scalar(
            conn, "SELECT COUNT(*) FROM drug_term_final WHERE mapping_stage = 'MANUAL_EXCLUDE'"
        )
        term_manual_mapped = scalar(
            conn,
            "SELECT COUNT(*) FROM drug_term_final WHERE COALESCE(manual_action,'') = 'MAP_TO_TERM' "
            "AND (COALESCE(rxcui,'') <> '' OR COALESCE(inchi,'') <> '')",
        )

        row_total = scalar(conn, "SELECT COALESCE(SUM(n_rows),0) FROM drug_term_final")
        row_mapped_rxcui = scalar(conn, "SELECT COALESCE(SUM(n_rows),0) FROM drug_term_final WHERE COALESCE(rxcui,'') <> ''")
        row_mapped_dbid = scalar(conn, "SELECT COALESCE(SUM(n_rows),0) FROM drug_term_final WHERE COALESCE(dbid,'') <> ''")
        row_mapped_inchi = scalar(conn, "SELECT COALESCE(SUM(n_rows),0) FROM drug_term_final WHERE COALESCE(inchi,'') <> ''")
        row_manual_excluded = scalar(
            conn, "SELECT COALESCE(SUM(n_rows),0) FROM drug_term_final WHERE mapping_stage = 'MANUAL_EXCLUDE'"
        )

        finished = dt.datetime.now().isoformat(timespec="seconds")
        summary = {
            "started_at": started,
            "finished_at": finished,
            "term_total": term_total,
            "term_mapped_rxcui": term_mapped_rxcui,
            "term_mapped_dbid": term_mapped_dbid,
            "term_mapped_inchi": term_mapped_inchi,
            "term_manual_excluded": term_manual_excluded,
            "term_manual_mapped": term_manual_mapped,
            "row_total": row_total,
            "row_mapped_rxcui": row_mapped_rxcui,
            "row_mapped_dbid": row_mapped_dbid,
            "row_mapped_inchi": row_mapped_inchi,
            "row_manual_excluded": row_manual_excluded,
            "term_coverage_rxcui": round(term_mapped_rxcui / term_total, 6) if term_total else 0.0,
            "term_coverage_dbid": round(term_mapped_dbid / term_total, 6) if term_total else 0.0,
            "term_coverage_inchi": round(term_mapped_inchi / term_total, 6) if term_total else 0.0,
            "row_coverage_rxcui": round(row_mapped_rxcui / row_total, 6) if row_total else 0.0,
            "row_coverage_dbid": round(row_mapped_dbid / row_total, 6) if row_total else 0.0,
            "row_coverage_inchi": round(row_mapped_inchi / row_total, 6) if row_total else 0.0,
            "rxnconso_scan": scan_stats,
            "rx_candidate_term_count": len(rx_by_term),
            "rxcui_with_drugbank_links": len(rxcui_to_dbids),
            "exported_final_rows": n_final,
        }

        (args.output_dir / "step3_report.json").write_text(
            json.dumps(summary, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )

        print("[DONE] Step3 completed", flush=True)
        print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    raise SystemExit(main())

