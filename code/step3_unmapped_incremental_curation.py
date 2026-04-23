#!/usr/bin/env python3
"""Step 3 supplementary curation for unmapped drug names.

Purpose:
- review high-frequency unmapped Step3 drug terms in fixed-size batches
- apply manual normalization choices before repeating RxNorm/DrugBank matching
- save incremental curation decisions for reproducibility

Inputs:
- outputs/step3/faers_step3.sqlite
- resources/step3_manual/

Outputs:
- outputs/step3/*.csv
- resources/step3_manual/*.csv
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
    r"\b(?:TABLET|TABLETS|TAB|CAP|CAPSULE|CAPSULES|INJECTION|INJECTABLE|SOLUTION|SUSPENSION|SYRUP|CREAM|OINTMENT|"
    r"PATCH|UNKNOWN|UNK|NOS|FORMULATION|GENERIC|BLINDED|DOSE|DOSES|MG|MCG|ML|IU|VACCINE|NDA)\b"
)
NON_DRUG_HINT_RE = re.compile(
    r"\b(?:UNKNOWN|UNSPECIFIED|NO DRUG NAME|NO CONCURRENT|THERAPY|RADIATION|RADIOTHERAPY|DEVICE|HOMEOPATH|"
    r"INVESTIGATIONAL|CHEMOTHERAPY|BIRTH CONTROL PILLS|DRUG USED IN DIABETES|BLOOD PRESSURE MEDICATION|"
    r"HUMAN RED BLOOD CELL|HUMAN PLATELET|PLATELETS)\b"
)

SALT_SUFFIX_RE = re.compile(
    r"\b(?:HYDROCHLORIDE|HCL|BESYLATE|MALEATE|MESYLATE|FUMARATE|PHOSPHATE|DIPHOSPHATE|"
    r"SODIUM|POTASSIUM|CALCIUM|MAGNESIUM|SUCCINATE|TARTRATE|CITRATE|ACETATE|NITRATE|"
    r"MONOHYDRATE|DIHYDRATE|TRIHYDRATE|ANHYDROUS)\b$"
)

ALLOWED_RELA = {
    "HAS_INGREDIENT",
    "INGREDIENT_OF",
    "HAS_PRECISE_INGREDIENT",
    "PRECISE_INGREDIENT_OF",
    "TRADENAME_OF",
    "HAS_TRADENAME",
    "CONSISTS_OF",
    "CONSTITUTES",
    "FORM_OF",
    "HAS_FORM",
}


@dataclass(frozen=True)
class RxCandidate:
    rxcui: str
    sab: str
    tty: str
    code: str
    string: str


@dataclass
class MapResult:
    success: bool
    stage: str
    target_term: str
    rxcui: str
    sab: str
    tty: str
    code: str
    rx_string: str
    dbid: str
    drugbank_name: str
    inchi: str
    candidate_count: int
    method: str


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


def strip_salt_suffix(s: str) -> str:
    t = normalize_drug_term(s)
    if not t:
        return ""
    prev = t
    while True:
        nxt = SALT_SUFFIX_RE.sub("", prev).strip()
        nxt = re.sub(r"\s+", " ", nxt).strip()
        if not nxt or nxt == prev:
            break
        prev = nxt
    return prev


def build_variants(term_std: str) -> List[Tuple[int, str, str]]:
    variants: List[Tuple[int, str, str]] = []
    seen: Set[str] = set()

    def add(pri: int, stage: str, value: str) -> None:
        v = normalize_drug_term(value)
        if not v or v in seen:
            return
        seen.add(v)
        variants.append((pri, stage, v))

    base = normalize_drug_term(term_std)
    add(1, "EXACT", base)
    add(2, "CLEANED", clean_drug_term(term_std))
    if " AND " in base:
        add(3, "AND_TO_SLASH", base.replace(" AND ", " / "))
    if " OR " in base:
        add(3, "OR_TO_SLASH", base.replace(" OR ", " / "))
    for m in re.finditer(r"\(([^()]{2,})\)", base):
        add(4, "PAREN", m.group(1))
    for sep in ["/", "+", ",", ";", "&"]:
        if sep in base:
            for part in base.split(sep):
                add(5, "SPLIT", part)
    return sorted(variants, key=lambda x: (x[0], x[1], x[2]))


def candidate_rank(c: RxCandidate, target_norm: str) -> Tuple[int, int, int, int, int]:
    sab_p = SAB_PRIORITY.get(c.sab, 99)
    tty_p = TTY_PRIORITY.get(c.tty, 99)
    s_norm = normalize_drug_term(c.string)
    exact = 0 if s_norm == target_norm else 1
    len_diff = abs(len(s_norm) - len(target_norm))
    return (sab_p, tty_p, exact, len_diff, len(s_norm))


def scan_rxnconso(
    rxnconso_path: Path,
) -> Tuple[Dict[str, List[RxCandidate]], Dict[str, List[RxCandidate]], Dict[str, Set[str]], Dict[str, int]]:
    rx_exact: Dict[str, List[RxCandidate]] = defaultdict(list)
    rx_clean: Dict[str, List[RxCandidate]] = defaultdict(list)
    rxcui_to_dbids: Dict[str, Set[str]] = defaultdict(set)
    stats = {"lines_total": 0, "lines_english_active": 0, "drugbank_links": 0}

    with rxnconso_path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        for line in f:
            stats["lines_total"] += 1
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
                stats["drugbank_links"] += 1

            if lat != "ENG" or suppress != "N":
                continue
            stats["lines_english_active"] += 1
            cand = RxCandidate(rxcui=rxcui, sab=sab, tty=tty, code=code, string=string)
            ne = normalize_drug_term(string)
            nc = clean_drug_term(string)
            if ne:
                rx_exact[ne].append(cand)
            if nc:
                rx_clean[nc].append(cand)

    for d in (rx_exact, rx_clean):
        for key, vals in list(d.items()):
            seen = set()
            uniq = []
            for c in vals:
                k = (c.rxcui, c.sab, c.tty, c.code, c.string)
                if k in seen:
                    continue
                seen.add(k)
                uniq.append(c)
            d[key] = uniq
    return rx_exact, rx_clean, rxcui_to_dbids, stats


def relation_allowed(rela: str, rel: str) -> bool:
    rla = normalize_text(rela, upper=True)
    if rla and rla in ALLOWED_RELA:
        return True
    r = normalize_text(rel, upper=True)
    if r in {"RN", "RB"}:
        return True
    return False


def collect_candidate_rxcuis(
    selected_terms: Sequence[Dict[str, str]],
    rx_exact: Dict[str, List[RxCandidate]],
    rx_clean: Dict[str, List[RxCandidate]],
) -> Set[str]:
    out: Set[str] = set()
    for row in selected_terms:
        term_std = row.get("term_std", "") or ""
        variants = build_variants(term_std)
        for _, _, norm_v in variants:
            for c in rx_exact.get(norm_v, []):
                if c.rxcui:
                    out.add(c.rxcui)
            ck = clean_drug_term(norm_v)
            if ck:
                for c in rx_clean.get(ck, []):
                    if c.rxcui:
                        out.add(c.rxcui)
    return out


def build_rxcui_bridge_map(
    rxnrel_path: Path,
    candidate_rxcuis: Set[str],
    rxcui_with_dbid: Set[str],
) -> Dict[str, Set[str]]:
    if not rxnrel_path.exists() or not candidate_rxcuis or not rxcui_with_dbid:
        return {}

    direct_bridge: Dict[str, Set[str]] = defaultdict(set)
    neighbors: Dict[str, Set[str]] = defaultdict(set)
    neighbor_set: Set[str] = set()

    # Pass 1: candidate -> (dbid_rxcui or neighbor)
    with rxnrel_path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        for line in f:
            parts = line.rstrip("\n\r").split("|")
            if len(parts) < 15:
                continue
            r1 = parts[0].strip()
            rel = parts[3].strip()
            r2 = parts[4].strip()
            rela = parts[7].strip()
            suppress = parts[14].strip()
            if suppress and suppress != "N":
                continue
            if not relation_allowed(rela, rel):
                continue

            if r1 in candidate_rxcuis:
                if r2 in rxcui_with_dbid:
                    direct_bridge[r1].add(r2)
                else:
                    neighbors[r1].add(r2)
                    neighbor_set.add(r2)
            if r2 in candidate_rxcuis:
                if r1 in rxcui_with_dbid:
                    direct_bridge[r2].add(r1)
                else:
                    neighbors[r2].add(r1)
                    neighbor_set.add(r1)

    if not neighbor_set:
        return {k: set(v) for k, v in direct_bridge.items()}

    neighbor_to_dbid: Dict[str, Set[str]] = defaultdict(set)

    # Pass 2: neighbor -> dbid_rxcui
    with rxnrel_path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        for line in f:
            parts = line.rstrip("\n\r").split("|")
            if len(parts) < 15:
                continue
            r1 = parts[0].strip()
            rel = parts[3].strip()
            r2 = parts[4].strip()
            rela = parts[7].strip()
            suppress = parts[14].strip()
            if suppress and suppress != "N":
                continue
            if not relation_allowed(rela, rel):
                continue

            if r1 in neighbor_set and r2 in rxcui_with_dbid:
                neighbor_to_dbid[r1].add(r2)
            if r2 in neighbor_set and r1 in rxcui_with_dbid:
                neighbor_to_dbid[r2].add(r1)

    out: Dict[str, Set[str]] = defaultdict(set)
    for c in candidate_rxcuis:
        out[c].update(direct_bridge.get(c, set()))
        for n in neighbors.get(c, set()):
            out[c].update(neighbor_to_dbid.get(n, set()))
    return {k: v for k, v in out.items() if v}


def load_drugbank_maps(drugbank_dir: Path) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, Set[str]], Dict[str, Set[str]]]:
    dname_path = drugbank_dir / "XML_dbid_dname.csv"
    inchi_path = drugbank_dir / "XML_dbid_inchi.csv"
    dbname: Dict[str, str] = {}
    dbinchi: Dict[str, str] = {}
    db_exact: Dict[str, Set[str]] = defaultdict(set)
    db_clean: Dict[str, Set[str]] = defaultdict(set)

    with dname_path.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            dbid = normalize_text(row.get("DrugID", ""), upper=True)
            name_raw = row.get("DrugName", "") or ""
            if not dbid:
                continue
            dbname.setdefault(dbid, normalize_text(name_raw, upper=False))
            ne = normalize_drug_term(name_raw)
            nc = clean_drug_term(name_raw)
            if ne:
                db_exact[ne].add(dbid)
            if nc:
                db_clean[nc].add(dbid)

    with inchi_path.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            dbid = normalize_text(row.get("DrugID", ""), upper=True)
            inchi = normalize_text(row.get("InChI", ""), upper=False)
            if dbid and inchi:
                dbinchi.setdefault(dbid, inchi)
    return dbname, dbinchi, db_exact, db_clean


def choose_dbid_for_rxcui(
    rxcui: str,
    rxcui_to_dbids: Dict[str, Set[str]],
    rxcui_bridge_to_dbid_rxcui: Dict[str, Set[str]],
    dbname: Dict[str, str],
    dbinchi: Dict[str, str],
) -> Tuple[str, str, str]:
    dbids = sorted(rxcui_to_dbids.get(rxcui, set()))
    if not dbids:
        for related in sorted(rxcui_bridge_to_dbid_rxcui.get(rxcui, set())):
            dbids.extend(sorted(rxcui_to_dbids.get(related, set())))
        if dbids:
            dbids = sorted(set(dbids))
    if not dbids:
        return "", "", ""
    with_inchi = [d for d in dbids if d in dbinchi]
    chosen = with_inchi[0] if with_inchi else dbids[0]
    return chosen, dbname.get(chosen, ""), dbinchi.get(chosen, "")


def choose_dbid_by_name(
    terms: Iterable[str],
    db_exact: Dict[str, Set[str]],
    db_clean: Dict[str, Set[str]],
    dbname: Dict[str, str],
    dbinchi: Dict[str, str],
) -> Tuple[str, str, str]:
    for t in terms:
        ne = normalize_drug_term(t)
        nc = clean_drug_term(t)
        for key, d in ((ne, db_exact), (nc, db_clean)):
            if not key:
                continue
            for dbid in sorted(d.get(key, set())):
                if dbid in dbinchi:
                    return dbid, dbname.get(dbid, ""), dbinchi.get(dbid, "")
    return "", "", ""


def try_map_term(
    term: str,
    rx_exact: Dict[str, List[RxCandidate]],
    rx_clean: Dict[str, List[RxCandidate]],
    rxcui_to_dbids: Dict[str, Set[str]],
    rxcui_bridge_to_dbid_rxcui: Dict[str, Set[str]],
    db_exact: Dict[str, Set[str]],
    db_clean: Dict[str, Set[str]],
    dbname: Dict[str, str],
    dbinchi: Dict[str, str],
) -> MapResult:
    variants = build_variants(term)
    for _, stage, norm_v in variants:
        pool: List[RxCandidate] = []
        pool.extend(rx_exact.get(norm_v, []))
        ckey = clean_drug_term(norm_v)
        if ckey:
            pool.extend(rx_clean.get(ckey, []))
        if not pool:
            continue
        seen = set()
        uniq = []
        for c in pool:
            k = (c.rxcui, c.sab, c.tty, c.code, c.string)
            if k in seen:
                continue
            seen.add(k)
            uniq.append(c)
        best = sorted(uniq, key=lambda c: candidate_rank(c, normalize_drug_term(term)))[0]
        dbid, dbnm, inchi = choose_dbid_for_rxcui(
            best.rxcui,
            rxcui_to_dbids,
            rxcui_bridge_to_dbid_rxcui,
            dbname,
            dbinchi,
        )
        if not inchi:
            dbid, dbnm, inchi = choose_dbid_by_name(
                [term, best.string, clean_drug_term(term), clean_drug_term(best.string), strip_salt_suffix(term)],
                db_exact,
                db_clean,
                dbname,
                dbinchi,
            )
        if inchi or best.rxcui:
            return MapResult(
                success=True,
                stage=stage,
                target_term=term,
                rxcui=best.rxcui,
                sab=best.sab,
                tty=best.tty,
                code=best.code,
                rx_string=best.string,
                dbid=dbid,
                drugbank_name=dbnm,
                inchi=inchi,
                candidate_count=len(uniq),
                method="RXNORM",
            )
    dbid, dbnm, inchi = choose_dbid_by_name([term], db_exact, db_clean, dbname, dbinchi)
    if inchi:
        return MapResult(
            success=True,
            stage="DB_NAME_EXACT",
            target_term=term,
            rxcui="",
            sab="",
            tty="",
            code="",
            rx_string="",
            dbid=dbid,
            drugbank_name=dbnm,
            inchi=inchi,
            candidate_count=1,
            method="DRUGBANK_NAME",
        )
    return MapResult(False, "UNMAPPED", term, "", "", "", "", "", "", "", "", 0, "NONE")


def is_non_drug_term(term: str) -> bool:
    t = normalize_drug_term(term)
    if not t:
        return True
    if t in {".", "-", "?", "??", "N/A"}:
        return True
    return bool(NON_DRUG_HINT_RE.search(t))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Incremental curation for step3 unmapped terms")
    project_root = Path(__file__).resolve().parents[1]
    data_root = project_root.parents[1]
    p.add_argument("--base-final-csv", type=Path, default=project_root / "outputs" / "step3" / "drug_term_final.csv")
    p.add_argument("--unmapped-csv", type=Path, default=project_root / "outputs" / "step3" / "unmapped_drug_terms_top.csv")
    p.add_argument("--rxnconso", type=Path, default=data_root / "reference_data" / "RxNorm_full_03022026" / "rrf" / "RXNCONSO.RRF")
    p.add_argument("--rxnrel", type=Path, default=data_root / "reference_data" / "RxNorm_full_03022026" / "rrf" / "RXNREL.RRF")
    p.add_argument("--drugbank-dir", type=Path, default=data_root / "reference_data" / "drugbank_5.1.15")
    p.add_argument("--manual-dir", type=Path, default=project_root / "resources" / "step3_manual")
    p.add_argument("--output-dir", type=Path, default=project_root / "outputs" / "step3")
    p.add_argument("--batch-size", type=int, default=300)
    p.add_argument("--top-n", type=int, default=5000)
    p.add_argument("--all-unmapped", action="store_true")
    p.add_argument("--progress-interval", type=int, default=5000)
    return p.parse_args()


def load_corrections(path: Path) -> Dict[str, Dict[str, str]]:
    out: Dict[str, Dict[str, str]] = {}
    if not path.exists():
        return out
    with path.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            key = normalize_drug_term(row.get("term_std", ""))
            if not key:
                continue
            out[key] = {
                "term_std": row.get("term_std", ""),
                "action": normalize_text(row.get("action", ""), upper=True),
                "target_term": row.get("target_term", ""),
                "note": row.get("note", ""),
            }
    return out


def main() -> int:
    args = parse_args()
    started = dt.datetime.now().isoformat(timespec="seconds")
    args.manual_dir.mkdir(parents=True, exist_ok=True)
    (args.manual_dir / "batches").mkdir(parents=True, exist_ok=True)

    print("[INFO] loading baseline csv ...", flush=True)
    base_rows: Dict[str, Dict[str, str]] = {}
    with args.base_final_csv.open("r", encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            base_rows[row["term_std"]] = row
    print(f"[INFO] baseline rows={len(base_rows)}", flush=True)

    print("[INFO] selecting unmapped terms ...", flush=True)
    top_unmapped: List[Dict[str, str]] = []
    if args.all_unmapped:
        for row in base_rows.values():
            if (row.get("inchi", "") or "").strip():
                continue
            if (row.get("mapping_stage", "") or "") == "MANUAL_EXCLUDE":
                continue
            top_unmapped.append(
                {
                    "term_std": row.get("term_std", ""),
                    "term_norm": row.get("term_norm", ""),
                    "n_rows": row.get("n_rows", "0"),
                    "n_cases": row.get("n_cases", "0"),
                }
            )
        top_unmapped.sort(key=lambda r: (-int(r.get("n_rows", "0") or 0), r.get("term_std", "")))
        if args.top_n > 0:
            top_unmapped = top_unmapped[: args.top_n]
    else:
        with args.unmapped_csv.open("r", encoding="utf-8", newline="") as f:
            for i, row in enumerate(csv.DictReader(f), start=1):
                if args.top_n > 0 and i > args.top_n:
                    break
                top_unmapped.append(row)
    print(f"[INFO] selected unmapped={len(top_unmapped)}", flush=True)

    master_path = args.manual_dir / "corrections_master.csv"
    corrections = load_corrections(master_path)
    seed_path = args.manual_dir / "corrections_batch1_top300.csv"
    for k, v in load_corrections(seed_path).items():
        corrections.setdefault(k, v)
    print(f"[INFO] preload corrections={len(corrections)}", flush=True)

    print("[INFO] loading DrugBank ...", flush=True)
    dbname, dbinchi, db_exact, db_clean = load_drugbank_maps(args.drugbank_dir)
    print(f"[INFO] DrugBank names={len(dbname)} inchi={len(dbinchi)}", flush=True)

    print("[INFO] scanning RXNCONSO ...", flush=True)
    rx_exact, rx_clean, rxcui_to_dbids, rx_stats = scan_rxnconso(args.rxnconso)
    print(f"[INFO] RXNCONSO eng_active={rx_stats['lines_english_active']}", flush=True)

    print("[INFO] building RXNREL bridge map ...", flush=True)
    candidate_rxcuis = collect_candidate_rxcuis(top_unmapped, rx_exact, rx_clean)
    rxcui_bridge = build_rxcui_bridge_map(
        args.rxnrel,
        candidate_rxcuis=candidate_rxcuis,
        rxcui_with_dbid=set(rxcui_to_dbids.keys()),
    )
    print(
        f"[INFO] RXNREL bridge ready: candidate_rxcui={len(candidate_rxcuis)} bridged={len(rxcui_bridge)}",
        flush=True,
    )

    decisions: List[Dict[str, str]] = []
    total = len(top_unmapped)
    batches = (total + args.batch_size - 1) // args.batch_size

    for b in range(batches):
        lo = b * args.batch_size
        hi = min(total, (b + 1) * args.batch_size)
        batch = top_unmapped[lo:hi]
        b_map = 0
        b_excl = 0
        b_unres = 0
        out_batch = []

        for row in batch:
            term_std = row["term_std"]
            key = normalize_drug_term(term_std)
            corr = corrections.get(key, None)
            action = ""
            target = ""
            note = ""
            status = "UNRESOLVED"

            if corr:
                action = corr.get("action", "")
                target = corr.get("target_term", "")
                note = corr.get("note", "")
            else:
                if is_non_drug_term(term_std):
                    action = "EXCLUDE_NON_DRUG"
                    note = "auto_non_drug_rule"
                else:
                    candidates: List[str] = []
                    seen = set()

                    def add(v: str) -> None:
                        n = normalize_drug_term(v)
                        if n and n not in seen:
                            seen.add(n)
                            candidates.append(n)

                    add(term_std)
                    add(clean_drug_term(term_std))
                    n = normalize_drug_term(term_std)
                    add(n.replace(" AND ", " / "))
                    add(n.replace(" OR ", " / "))
                    add(strip_salt_suffix(term_std))
                    add(strip_salt_suffix(clean_drug_term(term_std)))
                    for cand in candidates:
                        mr = try_map_term(
                            cand,
                            rx_exact,
                            rx_clean,
                            rxcui_to_dbids,
                            rxcui_bridge,
                            db_exact,
                            db_clean,
                            dbname,
                            dbinchi,
                        )
                        if mr.success and mr.inchi:
                            action = "MAP_TO_TERM"
                            target = cand
                            note = f"auto_{mr.method}_{mr.stage}"
                            break

            if action == "MAP_TO_TERM":
                mr = try_map_term(
                    target,
                    rx_exact,
                    rx_clean,
                    rxcui_to_dbids,
                    rxcui_bridge,
                    db_exact,
                    db_clean,
                    dbname,
                    dbinchi,
                )
                if mr.success and mr.inchi:
                    status = "MAP_SUCCESS"
                    b_map += 1
                    corrections[key] = {"term_std": term_std, "action": action, "target_term": target, "note": note}
                else:
                    action = ""
                    target = ""
                    note = ""
                    status = "UNRESOLVED"
                    b_unres += 1
            elif action == "EXCLUDE_NON_DRUG":
                status = "EXCLUDE"
                b_excl += 1
                corrections[key] = {"term_std": term_std, "action": action, "target_term": "", "note": note}
            else:
                b_unres += 1

            rec = {
                "term_std": term_std,
                "term_norm": row.get("term_norm", ""),
                "n_rows": row.get("n_rows", "0"),
                "n_cases": row.get("n_cases", "0"),
                "action": action,
                "target_term": target,
                "status": status,
                "note": note,
            }
            decisions.append(rec)
            out_batch.append(rec)

        batch_path = args.manual_dir / "batches" / f"batch_{b + 1:03d}.csv"
        with batch_path.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(out_batch[0].keys()) if out_batch else ["term_std"])
            w.writeheader()
            w.writerows(out_batch)
        print(f"[PROGRESS] batch {b + 1}/{batches} | processed {hi}/{total} | map={b_map} exclude={b_excl} unresolved={b_unres}", flush=True)

    print("[INFO] applying updates to baseline ...", flush=True)
    upd = Counter()
    for d in decisions:
        term_std = d["term_std"]
        action = d["action"]
        status = d["status"]
        target = d["target_term"]
        row = base_rows.get(term_std)
        if not row:
            continue
        row.setdefault("manual_action", "")
        row.setdefault("manual_target", "")

        if action == "EXCLUDE_NON_DRUG":
            row["manual_action"] = "EXCLUDE_NON_DRUG"
            row["manual_target"] = ""
            row["mapping_stage"] = "MANUAL_EXCLUDE"
            row["candidate_count"] = "0"
            row["rxcui"] = ""
            row["sab"] = ""
            row["tty"] = ""
            row["code"] = ""
            row["rx_string"] = ""
            row["dbid"] = ""
            row["drugbank_name"] = ""
            row["inchi"] = ""
            upd["exclude"] += 1
        elif action == "MAP_TO_TERM" and status == "MAP_SUCCESS":
            mr = try_map_term(
                target,
                rx_exact,
                rx_clean,
                rxcui_to_dbids,
                rxcui_bridge,
                db_exact,
                db_clean,
                dbname,
                dbinchi,
            )
            if mr.success and mr.inchi:
                row["manual_action"] = "MAP_TO_TERM"
                row["manual_target"] = target
                row["mapping_stage"] = f"MANUAL_{mr.stage}"
                row["candidate_count"] = str(mr.candidate_count)
                row["rxcui"] = mr.rxcui
                row["sab"] = mr.sab
                row["tty"] = mr.tty
                row["code"] = mr.code
                row["rx_string"] = mr.rx_string
                row["dbid"] = mr.dbid
                row["drugbank_name"] = mr.drugbank_name
                row["inchi"] = mr.inchi
                upd["map_success"] += 1
            else:
                upd["map_failed"] += 1
        else:
            upd["unresolved"] += 1
    print(f"[INFO] update stats={dict(upd)}", flush=True)

    with master_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["term_std", "action", "target_term", "note"])
        w.writeheader()
        for _, v in sorted(corrections.items(), key=lambda kv: kv[1]["term_std"]):
            w.writerow(v)

    decisions_path = args.manual_dir / "decisions_top5000.csv"
    with decisions_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(decisions[0].keys()) if decisions else ["term_std"])
        w.writeheader()
        w.writerows(decisions)

    out_final = args.output_dir / "drug_term_final.csv"
    out_unmapped = args.output_dir / "unmapped_drug_terms_top.csv"
    out_stage = args.output_dir / "term_stage_counts.csv"
    out_report = args.output_dir / "step3_report.json"
    out_db = args.output_dir / "faers_step3.sqlite"
    backup_csv = args.output_dir / "drug_term_final.backup_before_incremental.csv"
    if not backup_csv.exists():
        backup_csv.write_bytes(args.base_final_csv.read_bytes())

    fields = ["term_std", "term_norm", "n_rows", "n_cases", "manual_action", "manual_target", "mapping_stage", "candidate_count", "rxcui", "sab", "tty", "code", "rx_string", "dbid", "drugbank_name", "inchi"]
    with out_final.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for _, row in sorted(base_rows.items(), key=lambda kv: (-int(kv[1].get("n_rows", "0") or 0), kv[0])):
            w.writerow({k: row.get(k, "") for k in fields})

    sorted_rows = sorted(base_rows.values(), key=lambda r: (-int(r.get("n_rows", "0") or 0), r.get("term_std", "")))
    output_unmapped_limit = args.top_n if args.top_n and args.top_n > 0 else 5000
    unmapped_rows = [
        r for r in sorted_rows if not (r.get("inchi", "") or "").strip() and r.get("mapping_stage", "") != "MANUAL_EXCLUDE"
    ][: output_unmapped_limit]
    unmapped_output_path = out_unmapped
    try:
        with out_unmapped.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["term_std", "term_norm", "n_rows", "n_cases"])
            w.writeheader()
            for r in unmapped_rows:
                w.writerow({"term_std": r.get("term_std", ""), "term_norm": r.get("term_norm", ""), "n_rows": r.get("n_rows", "0"), "n_cases": r.get("n_cases", "0")})
    except PermissionError:
        fallback = args.output_dir / "unmapped_drug_terms_top.locked_fallback.csv"
        with fallback.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["term_std", "term_norm", "n_rows", "n_cases"])
            w.writeheader()
            for r in unmapped_rows:
                w.writerow({"term_std": r.get("term_std", ""), "term_norm": r.get("term_norm", ""), "n_rows": r.get("n_rows", "0"), "n_cases": r.get("n_cases", "0")})
        unmapped_output_path = fallback
        print(f"[WARN] target unmapped csv is locked, wrote fallback: {fallback}", flush=True)

    st_term = Counter()
    st_row = Counter()
    for r in base_rows.values():
        s = r.get("mapping_stage", "") or "UNMAPPED"
        n = int(r.get("n_rows", "0") or 0)
        st_term[s] += 1
        st_row[s] += n
    with out_stage.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mapping_stage", "term_count", "row_count"])
        for s in sorted(st_term.keys(), key=lambda x: (-st_row[x], x)):
            w.writerow([s, st_term[s], st_row[s]])

    if out_db.exists():
        out_db.unlink()
    conn = sqlite3.connect(out_db)
    try:
        conn.execute("PRAGMA journal_mode = WAL;")
        conn.execute("PRAGMA synchronous = NORMAL;")
        conn.execute(
            """
            CREATE TABLE drug_term_final (
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
            )
            """
        )
        conn.execute("CREATE INDEX idx_dtf_inchi ON drug_term_final(inchi)")
        ins = "INSERT INTO drug_term_final(term_std,term_norm,n_rows,n_cases,manual_action,manual_target,mapping_stage,candidate_count,rxcui,sab,tty,code,rx_string,dbid,drugbank_name,inchi) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
        rows = []
        for r in base_rows.values():
            rows.append(
                (
                    r.get("term_std", ""),
                    r.get("term_norm", ""),
                    int(r.get("n_rows", "0") or 0),
                    int(r.get("n_cases", "0") or 0),
                    r.get("manual_action", ""),
                    r.get("manual_target", ""),
                    r.get("mapping_stage", ""),
                    int(r.get("candidate_count", "0") or 0),
                    r.get("rxcui", ""),
                    r.get("sab", ""),
                    r.get("tty", ""),
                    r.get("code", ""),
                    r.get("rx_string", ""),
                    r.get("dbid", ""),
                    r.get("drugbank_name", ""),
                    r.get("inchi", ""),
                )
            )
        conn.executemany(ins, rows)
        conn.commit()
    finally:
        conn.close()

    term_total = len(base_rows)
    row_total = sum(int(r.get("n_rows", "0") or 0) for r in base_rows.values())
    term_mapped_rxcui = sum(1 for r in base_rows.values() if (r.get("rxcui", "") or "").strip())
    term_mapped_dbid = sum(1 for r in base_rows.values() if (r.get("dbid", "") or "").strip())
    term_mapped_inchi = sum(1 for r in base_rows.values() if (r.get("inchi", "") or "").strip())
    row_mapped_rxcui = sum(int(r.get("n_rows", "0") or 0) for r in base_rows.values() if (r.get("rxcui", "") or "").strip())
    row_mapped_dbid = sum(int(r.get("n_rows", "0") or 0) for r in base_rows.values() if (r.get("dbid", "") or "").strip())
    row_mapped_inchi = sum(int(r.get("n_rows", "0") or 0) for r in base_rows.values() if (r.get("inchi", "") or "").strip())
    term_manual_excluded = sum(1 for r in base_rows.values() if r.get("mapping_stage", "") == "MANUAL_EXCLUDE")
    row_manual_excluded = sum(int(r.get("n_rows", "0") or 0) for r in base_rows.values() if r.get("mapping_stage", "") == "MANUAL_EXCLUDE")

    mode_name = "incremental_unmapped_all" if args.all_unmapped else "incremental_unmapped_topN"
    report = {
        "started_at": started,
        "finished_at": dt.datetime.now().isoformat(timespec="seconds"),
        "mode": mode_name,
        "batch_size": args.batch_size,
        "batches": batches,
        "selected_unmapped_terms": total,
        "update_stats": dict(upd),
        "term_total": term_total,
        "term_mapped_rxcui": term_mapped_rxcui,
        "term_mapped_dbid": term_mapped_dbid,
        "term_mapped_inchi": term_mapped_inchi,
        "term_manual_excluded": term_manual_excluded,
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
        "rxnconso_scan": rx_stats,
        "unmapped_output_path": str(unmapped_output_path),
    }
    out_report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print("[DONE] incremental curation finished", flush=True)
    print(json.dumps(report, ensure_ascii=False, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

