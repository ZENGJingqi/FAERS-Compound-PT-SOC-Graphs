from __future__ import annotations

import json
import sqlite3
import time
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "outputs" / "step6_core_graphs" / "core_graph_full"
OUT_DB = OUT_DIR / "step6_compound_pt_soc_core_full.sqlite"
STEP2_DB = ROOT / "outputs" / "step2" / "faers_step2.sqlite"
STEP4_DB = ROOT / "outputs" / "step4" / "faers_step4.sqlite"

FETCH_SIZE = 100_000


def load_term_map() -> dict[str, list[str]]:
    conn = sqlite3.connect(STEP4_DB)
    cur = conn.cursor()
    term_map: dict[str, list[str]] = defaultdict(list)
    for term_norm, inchikey in cur.execute(
        """
        SELECT term_norm, inchikey
        FROM term_inchikey_map
        WHERE TRIM(COALESCE(inchikey, '')) <> ''
        """
    ):
        term_map[term_norm].append(inchikey)
    conn.close()
    return dict(term_map)


def init_role_table(conn: sqlite3.Connection, table: str) -> None:
    conn.execute(f"DROP TABLE IF EXISTS {table}")
    conn.execute(
        f"""
        CREATE TABLE {table} (
            primaryid TEXT NOT NULL,
            inchikey TEXT NOT NULL,
            PRIMARY KEY(primaryid, inchikey)
        ) WITHOUT ROWID
        """
    )
    conn.execute(f"CREATE INDEX idx_{table}_inchikey ON {table}(inchikey)")
    conn.commit()


def scalar(conn: sqlite3.Connection, q: str):
    return conn.execute(q).fetchone()[0]


def materialize_role(
    src_conn: sqlite3.Connection,
    out_conn: sqlite3.Connection,
    term_map: dict[str, list[str]],
    role: str,
) -> dict:
    table = f"report_compound_{role.lower()}"
    init_role_table(out_conn, table)

    total_rows = scalar(
        src_conn,
        f"SELECT COUNT(*) FROM case_drug WHERE role_cod='{role}' AND TRIM(COALESCE(drug_key_std,''))<>''",
    )
    cur = src_conn.cursor()
    cur.execute(
        """
        SELECT primaryid, drug_key_std
        FROM case_drug
        WHERE role_cod = ?
          AND TRIM(COALESCE(drug_key_std, '')) <> ''
        """,
        (role,),
    )

    processed = 0
    matched_source_rows = 0
    t0 = time.time()
    batch_no = 0
    while True:
        rows = cur.fetchmany(FETCH_SIZE)
        if not rows:
            break
        batch_no += 1
        processed += len(rows)
        pair_set: set[tuple[str, str]] = set()
        for primaryid, drug_key_std in rows:
            keys = term_map.get(drug_key_std)
            if not keys:
                continue
            matched_source_rows += 1
            for inchikey in keys:
                pair_set.add((primaryid, inchikey))
        if pair_set:
            out_conn.executemany(
                f"INSERT OR IGNORE INTO {table}(primaryid, inchikey) VALUES (?, ?)",
                list(pair_set),
            )
            out_conn.commit()
        if batch_no % 10 == 0:
            current = scalar(out_conn, f"SELECT COUNT(*) FROM {table}")
            print(
                f"[INFO] {role}: processed {processed}/{total_rows} source rows; unique report-compound rows={current}; elapsed={round(time.time()-t0,1)}s",
                flush=True,
            )

    return {
        "role": role,
        "source_rows": int(total_rows),
        "matched_source_rows": int(matched_source_rows),
        "unique_rows": int(scalar(out_conn, f"SELECT COUNT(*) FROM {table}")),
        "unique_reports": int(scalar(out_conn, f"SELECT COUNT(DISTINCT primaryid) FROM {table}")),
        "unique_compounds": int(scalar(out_conn, f"SELECT COUNT(DISTINCT inchikey) FROM {table}")),
        "elapsed_sec": round(time.time() - t0, 1),
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if OUT_DB.exists():
        OUT_DB.unlink()

    print("[INFO] loading term_norm -> inchikey map ...", flush=True)
    term_map = load_term_map()
    print(f"[INFO] loaded mapped terms: {len(term_map)}", flush=True)

    src_conn = sqlite3.connect(STEP2_DB)
    out_conn = sqlite3.connect(OUT_DB)
    out_conn.execute("PRAGMA journal_mode=WAL;")
    out_conn.execute("PRAGMA synchronous=OFF;")
    out_conn.execute("PRAGMA temp_store=MEMORY;")

    ps = materialize_role(src_conn, out_conn, term_map, "PS")
    ss = materialize_role(src_conn, out_conn, term_map, "SS")

    summary = {"PS": ps, "SS": ss}
    (OUT_DIR / "report_compound_role_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    src_conn.close()
    out_conn.close()


if __name__ == "__main__":
    main()
