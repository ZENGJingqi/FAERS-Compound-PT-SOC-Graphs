#!/usr/bin/env python3
"""Step 1: FAERS/AERS case-level ingestion, deduplication, and delete-list filtering.

This script builds a durable SQLite-backed case index from quarterly ZIP files.
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
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

ZIP_RE = re.compile(r"(?i)^(?:f?aers)_ascii_(20\d{2})[qQ]([1-4])\.zip$")
DEMO_ENTRY_RE = re.compile(r"(?i)(^|/)(demo\d{2}q[1-4](?:_new)?\.txt)$")
DELETE_ENTRY_RE = re.compile(r"(?i)(^|/)(delete\d{2}q[1-4]\.txt)$")

FIELD_SYNONYMS = {
    "primaryid": ("primaryid", "isr"),
    "caseid": ("caseid", "case"),
    "caseversion": ("caseversion", "foll_seq"),
    "i_f_code": ("i_f_code", "i_f_cod"),
    "fda_dt": ("fda_dt",),
    "event_dt": ("event_dt",),
    "sex": ("sex", "gndr_cod"),
    "age": ("age",),
    "reporter_country": ("reporter_country",),
}


@dataclass(frozen=True)
class ZipInfo:
    path: Path
    year: int
    quarter: int


@dataclass
class ZipStats:
    zip_name: str
    year: int
    quarter: int
    demo_entry: str
    delete_entry: str
    demo_rows_total: int = 0
    demo_rows_kept: int = 0
    skipped_missing_ids: int = 0
    malformed_rows_fixed: int = 0
    deleted_ids_added: int = 0


def normalize(s: Optional[str]) -> str:
    if s is None:
        return ""
    return s.strip()


def collect_zip_files(data_root: Path) -> List[ZipInfo]:
    infos: List[ZipInfo] = []
    for path in data_root.glob("*aers_ascii_*.zip"):
        match = ZIP_RE.match(path.name)
        if not match:
            continue
        infos.append(ZipInfo(path=path, year=int(match.group(1)), quarter=int(match.group(2))))
    infos.sort(key=lambda x: (x.year, x.quarter, x.path.name.lower()))
    return infos


def find_entry(archive: zipfile.ZipFile, pattern: re.Pattern[str]) -> Optional[zipfile.ZipInfo]:
    for entry in archive.infolist():
        if pattern.search(entry.filename):
            return entry
    return None


def build_header_index(header: Sequence[str]) -> Dict[str, Optional[int]]:
    lowered = [normalize(col).lower() for col in header]
    idx_map: Dict[str, Optional[int]] = {}
    for canonical, aliases in FIELD_SYNONYMS.items():
        idx = None
        for alias in aliases:
            if alias in lowered:
                idx = lowered.index(alias)
                break
        idx_map[canonical] = idx
    return idx_map


def cell(row: Sequence[str], idx: Optional[int]) -> str:
    if idx is None:
        return ""
    if idx < 0 or idx >= len(row):
        return ""
    return normalize(row[idx])


def parse_delete_ids(archive: zipfile.ZipFile, entry: Optional[zipfile.ZipInfo]) -> List[str]:
    if entry is None:
        return []
    out: List[str] = []
    with archive.open(entry) as raw:
        text = io.TextIOWrapper(raw, encoding="latin-1", errors="replace", newline="")
        for line in text:
            token = normalize(line)
            if not token:
                continue
            if "$" in token:
                token = normalize(token.split("$", 1)[0])
            if token.isdigit():
                out.append(token)
    return out


def iter_demo_rows(
    archive: zipfile.ZipFile,
    entry: zipfile.ZipInfo,
    year: int,
    quarter: int,
    zip_name: str,
    stats: ZipStats,
) -> Iterable[Tuple[str, str, str, str, str, str, int, int, str, str, str, str]]:
    with archive.open(entry) as raw:
        text = io.TextIOWrapper(raw, encoding="latin-1", errors="replace", newline="")
        reader = csv.reader(text, delimiter="$", quotechar='"')

        try:
            header = next(reader)
        except StopIteration:
            return

        idx_map = build_header_index(header)

        if idx_map["primaryid"] is None and idx_map["caseid"] is None:
            raise RuntimeError(f"No primary identifier columns found in {zip_name}:{entry.filename}")

        header_len = len(header)
        for row in reader:
            stats.demo_rows_total += 1
            if not row:
                stats.skipped_missing_ids += 1
                continue

            if len(row) < header_len:
                row = list(row) + [""] * (header_len - len(row))
                stats.malformed_rows_fixed += 1
            elif len(row) > header_len:
                row = list(row[: header_len - 1]) + ["$".join(row[header_len - 1 :])]
                stats.malformed_rows_fixed += 1

            primaryid = cell(row, idx_map["primaryid"])
            caseid = cell(row, idx_map["caseid"])
            caseversion = cell(row, idx_map["caseversion"])
            i_f_code = cell(row, idx_map["i_f_code"])
            fda_dt = cell(row, idx_map["fda_dt"])
            event_dt = cell(row, idx_map["event_dt"])
            sex = cell(row, idx_map["sex"])
            age = cell(row, idx_map["age"])
            reporter_country = cell(row, idx_map["reporter_country"])

            if not primaryid and not caseid:
                stats.skipped_missing_ids += 1
                continue
            if not caseid:
                caseid = primaryid
            if not primaryid:
                primaryid = caseid

            stats.demo_rows_kept += 1
            yield (
                caseid,
                primaryid,
                caseversion,
                i_f_code,
                fda_dt,
                event_dt,
                year,
                quarter,
                zip_name,
                sex,
                age,
                reporter_country,
            )


def init_db(conn: sqlite3.Connection, rebuild: bool) -> None:
    cur = conn.cursor()
    conn.execute("PRAGMA journal_mode = WAL;")
    conn.execute("PRAGMA synchronous = NORMAL;")
    conn.execute("PRAGMA temp_store = MEMORY;")
    conn.execute("PRAGMA cache_size = -200000;")

    if rebuild:
        cur.executescript(
            """
            DROP TABLE IF EXISTS demo_raw;
            DROP TABLE IF EXISTS deleted_ids;
            DROP TABLE IF EXISTS demo_scored;
            DROP TABLE IF EXISTS demo_primary_dedup;
            DROP TABLE IF EXISTS case_selected;
            """
        )

    cur.executescript(
        """
        CREATE TABLE IF NOT EXISTS demo_raw (
            caseid TEXT NOT NULL,
            primaryid TEXT NOT NULL,
            caseversion TEXT,
            i_f_code TEXT,
            fda_dt TEXT,
            event_dt TEXT,
            year INTEGER NOT NULL,
            quarter INTEGER NOT NULL,
            zip_name TEXT NOT NULL,
            sex TEXT,
            age TEXT,
            reporter_country TEXT
        );

        CREATE TABLE IF NOT EXISTS deleted_ids (
            id TEXT PRIMARY KEY,
            first_seen_zip TEXT
        );
        """
    )
    conn.commit()


def ingest_all(
    conn: sqlite3.Connection,
    zip_infos: Sequence[ZipInfo],
    batch_size: int = 50000,
) -> List[ZipStats]:
    stats_rows: List[ZipStats] = []
    insert_sql = (
        "INSERT INTO demo_raw(caseid,primaryid,caseversion,i_f_code,fda_dt,event_dt,year,quarter,zip_name,sex,age,reporter_country) "
        "VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"
    )
    insert_delete_sql = "INSERT OR IGNORE INTO deleted_ids(id, first_seen_zip) VALUES (?, ?)"

    for info in zip_infos:
        zip_name = info.path.name
        with zipfile.ZipFile(info.path, "r") as archive:
            demo_entry = find_entry(archive, DEMO_ENTRY_RE)
            delete_entry = find_entry(archive, DELETE_ENTRY_RE)

            stat = ZipStats(
                zip_name=zip_name,
                year=info.year,
                quarter=info.quarter,
                demo_entry=demo_entry.filename if demo_entry else "",
                delete_entry=delete_entry.filename if delete_entry else "",
            )

            if demo_entry is None:
                stats_rows.append(stat)
                print(f"[WARN] DEMO file not found in {zip_name}", flush=True)
                continue

            deleted_ids = parse_delete_ids(archive, delete_entry)
            if deleted_ids:
                conn.executemany(insert_delete_sql, ((x, zip_name) for x in deleted_ids))
                stat.deleted_ids_added = len(deleted_ids)

            batch: List[Tuple[str, str, str, str, str, str, int, int, str, str, str, str]] = []
            for row in iter_demo_rows(
                archive=archive,
                entry=demo_entry,
                year=info.year,
                quarter=info.quarter,
                zip_name=zip_name,
                stats=stat,
            ):
                batch.append(row)
                if len(batch) >= batch_size:
                    conn.executemany(insert_sql, batch)
                    conn.commit()
                    batch.clear()

            if batch:
                conn.executemany(insert_sql, batch)
                conn.commit()

            stats_rows.append(stat)
            print(
                f"[INGEST] {zip_name}: demo_rows={stat.demo_rows_total}, kept={stat.demo_rows_kept}, "
                f"fixed={stat.malformed_rows_fixed}, deleted_ids={stat.deleted_ids_added}",
                flush=True,
            )

    return stats_rows


def build_dedup_tables(conn: sqlite3.Connection) -> None:
    cur = conn.cursor()
    cur.executescript(
        """
        DROP TABLE IF EXISTS demo_scored;
        DROP TABLE IF EXISTS demo_primary_dedup;
        DROP TABLE IF EXISTS case_selected;

        CREATE TABLE demo_scored AS
        SELECT
            caseid,
            primaryid,
            COALESCE(caseversion, '') AS caseversion,
            COALESCE(i_f_code, '') AS i_f_code,
            COALESCE(fda_dt, '') AS fda_dt,
            COALESCE(event_dt, '') AS event_dt,
            year,
            quarter,
            zip_name,
            COALESCE(sex, '') AS sex,
            COALESCE(age, '') AS age,
            COALESCE(reporter_country, '') AS reporter_country,
            (year * 10 + quarter) AS quarter_idx,
            CASE
                WHEN TRIM(caseversion) GLOB '[0-9][0-9]*' THEN CAST(TRIM(caseversion) AS INTEGER)
                ELSE 0
            END AS caseversion_n,
            CASE
                WHEN LENGTH(TRIM(fda_dt)) = 8
                    AND TRIM(fda_dt) GLOB '[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]'
                THEN CAST(TRIM(fda_dt) AS INTEGER)
                ELSE 0
            END AS fda_dt_n,
            CASE
                WHEN UPPER(TRIM(i_f_code)) = 'F' THEN 2
                WHEN UPPER(TRIM(i_f_code)) = 'I' THEN 1
                ELSE 0
            END AS if_rank,
            CASE
                WHEN LENGTH(TRIM(event_dt)) = 8
                    AND TRIM(event_dt) GLOB '[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]'
                THEN CAST(TRIM(event_dt) AS INTEGER)
                ELSE 0
            END AS event_dt_n,
            CASE
                WHEN TRIM(primaryid) GLOB '[0-9][0-9]*' THEN CAST(TRIM(primaryid) AS INTEGER)
                ELSE 0
            END AS primaryid_n
        FROM demo_raw r
        WHERE TRIM(caseid) <> ''
          AND TRIM(primaryid) <> ''
          AND NOT EXISTS (
            SELECT 1
            FROM deleted_ids d
            WHERE d.id = TRIM(r.caseid) OR d.id = TRIM(r.primaryid)
          );

        CREATE INDEX idx_demo_scored_primaryid ON demo_scored(primaryid);
        CREATE INDEX idx_demo_scored_caseid ON demo_scored(caseid);

        CREATE TABLE demo_primary_dedup AS
        SELECT
            caseid,
            primaryid,
            caseversion,
            i_f_code,
            fda_dt,
            event_dt,
            year,
            quarter,
            zip_name,
            sex,
            age,
            reporter_country,
            quarter_idx,
            caseversion_n,
            fda_dt_n,
            if_rank,
            event_dt_n,
            primaryid_n
        FROM (
            SELECT
                s.caseid,
                s.primaryid,
                s.caseversion,
                s.i_f_code,
                s.fda_dt,
                s.event_dt,
                s.year,
                s.quarter,
                s.zip_name,
                s.sex,
                s.age,
                s.reporter_country,
                s.quarter_idx,
                s.caseversion_n,
                s.fda_dt_n,
                s.if_rank,
                s.event_dt_n,
                s.primaryid_n,
                ROW_NUMBER() OVER (
                    PARTITION BY s.primaryid
                    ORDER BY s.quarter_idx DESC, s.fda_dt_n DESC, s.if_rank DESC,
                             s.event_dt_n DESC, s.caseversion_n DESC, s.primaryid_n DESC
                ) AS rn_primary
            FROM demo_scored s
        ) x
        WHERE x.rn_primary = 1;

        CREATE INDEX idx_demo_primary_dedup_caseid ON demo_primary_dedup(caseid);

        CREATE TABLE case_selected AS
        SELECT
            caseid,
            primaryid,
            caseversion,
            i_f_code,
            fda_dt,
            event_dt,
            year,
            quarter,
            zip_name,
            sex,
            age,
            reporter_country
        FROM (
            SELECT
                p.caseid,
                p.primaryid,
                p.caseversion,
                p.i_f_code,
                p.fda_dt,
                p.event_dt,
                p.year,
                p.quarter,
                p.zip_name,
                p.sex,
                p.age,
                p.reporter_country,
                ROW_NUMBER() OVER (
                    PARTITION BY p.caseid
                    ORDER BY p.caseversion_n DESC, p.fda_dt_n DESC, p.if_rank DESC,
                             p.event_dt_n DESC, p.quarter_idx DESC, p.primaryid_n DESC
                ) AS rn_case
            FROM demo_primary_dedup p
        ) y
        WHERE y.rn_case = 1;

        CREATE INDEX idx_case_selected_primaryid ON case_selected(primaryid);
        CREATE INDEX idx_case_selected_caseid ON case_selected(caseid);
        """
    )
    conn.commit()


def export_csv(conn: sqlite3.Connection, output_csv: Path) -> int:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    query = (
        "SELECT caseid, primaryid, caseversion, i_f_code, fda_dt, event_dt, "
        "year, quarter, zip_name, sex, age, reporter_country "
        "FROM case_selected ORDER BY year, quarter, caseid"
    )
    rows = 0
    with output_csv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "caseid",
                "primaryid",
                "caseversion",
                "i_f_code",
                "fda_dt",
                "event_dt",
                "year",
                "quarter",
                "zip_name",
                "sex",
                "age",
                "reporter_country",
            ]
        )
        for rec in conn.execute(query):
            writer.writerow(rec)
            rows += 1
    return rows


def write_zip_stats_csv(stats_rows: Sequence[ZipStats], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "zip_name",
                "year",
                "quarter",
                "demo_entry",
                "delete_entry",
                "demo_rows_total",
                "demo_rows_kept",
                "skipped_missing_ids",
                "malformed_rows_fixed",
                "deleted_ids_added",
            ]
        )
        for s in stats_rows:
            writer.writerow(
                [
                    s.zip_name,
                    s.year,
                    s.quarter,
                    s.demo_entry,
                    s.delete_entry,
                    s.demo_rows_total,
                    s.demo_rows_kept,
                    s.skipped_missing_ids,
                    s.malformed_rows_fixed,
                    s.deleted_ids_added,
                ]
            )


def scalar(conn: sqlite3.Connection, sql: str) -> int:
    return int(conn.execute(sql).fetchone()[0])


def build_summary(
    conn: sqlite3.Connection,
    zip_infos: Sequence[ZipInfo],
    stats_rows: Sequence[ZipStats],
    latest_case_count: int,
    started_at: str,
    finished_at: str,
) -> Dict[str, object]:
    return {
        "started_at": started_at,
        "finished_at": finished_at,
        "zip_files_total": len(zip_infos),
        "zip_files_with_demo": sum(1 for s in stats_rows if s.demo_entry),
        "zip_files_with_delete": sum(1 for s in stats_rows if s.delete_entry),
        "delete_ids_unique": scalar(conn, "SELECT COUNT(*) FROM deleted_ids"),
        "demo_rows_raw": scalar(conn, "SELECT COUNT(*) FROM demo_raw"),
        "demo_rows_after_delete_filter": scalar(conn, "SELECT COUNT(*) FROM demo_scored"),
        "rows_after_primaryid_dedup": scalar(conn, "SELECT COUNT(*) FROM demo_primary_dedup"),
        "rows_after_caseid_dedup": scalar(conn, "SELECT COUNT(*) FROM case_selected"),
        "latest_case_rows_exported": latest_case_count,
        "distinct_caseid_raw": scalar(conn, "SELECT COUNT(DISTINCT caseid) FROM demo_raw"),
        "distinct_primaryid_raw": scalar(conn, "SELECT COUNT(DISTINCT primaryid) FROM demo_raw"),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Step 1 FAERS/AERS case-level dedup pipeline")
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(r"D:\博士文件\TCMMKG\data\AEMS_FDA不良反应数据"),
        help="Directory containing quarterly *aers_ascii_*.zip files",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(r"D:\博士文件\TCMMKG\data\AEMS_FDA不良反应数据\smiles_adr_project\outputs\step1"),
        help="Directory for step1 outputs",
    )
    parser.add_argument(
        "--rebuild",
        action="store_true",
        help="Drop and rebuild Step 1 tables in SQLite before ingestion",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    data_root: Path = args.data_root
    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    if not data_root.exists():
        print(f"[ERROR] data-root does not exist: {data_root}", file=sys.stderr)
        return 1

    started_at = dt.datetime.now().isoformat(timespec="seconds")

    zip_infos = collect_zip_files(data_root)
    if not zip_infos:
        print(f"[ERROR] No FAERS/AERS zip files found in {data_root}", file=sys.stderr)
        return 1

    sqlite_path = output_dir / "faers_step1.sqlite"
    latest_csv_path = output_dir / "faers_latest_case_keys.csv"
    per_zip_csv_path = output_dir / "step1_ingest_per_zip.csv"
    summary_json_path = output_dir / "step1_case_dedup_report.json"

    conn = sqlite3.connect(sqlite_path)
    try:
        init_db(conn, rebuild=args.rebuild)

        print(f"[INFO] Ingesting {len(zip_infos)} zip files...", flush=True)
        stats_rows = ingest_all(conn, zip_infos)

        print("[INFO] Building dedup tables...", flush=True)
        build_dedup_tables(conn)

        print(f"[INFO] Exporting latest case keys: {latest_csv_path}", flush=True)
        latest_case_count = export_csv(conn, latest_csv_path)

        write_zip_stats_csv(stats_rows, per_zip_csv_path)

        finished_at = dt.datetime.now().isoformat(timespec="seconds")
        summary = build_summary(
            conn=conn,
            zip_infos=zip_infos,
            stats_rows=stats_rows,
            latest_case_count=latest_case_count,
            started_at=started_at,
            finished_at=finished_at,
        )
        summary_json_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

        print("[DONE] Step 1 completed.", flush=True)
        print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
        return 0
    finally:
        conn.close()


if __name__ == "__main__":
    raise SystemExit(main())
