# Step 1 Execution Notes

Run date: 2026-04-11
Script: `scripts/step1_faers_case_dedup.py`

## Command
```bash
python scripts/step1_faers_case_dedup.py --rebuild
```

## Outputs
- `outputs/step1/faers_step1.sqlite`
- `outputs/step1/faers_latest_case_keys.csv`
- `outputs/step1/step1_ingest_per_zip.csv`
- `outputs/step1/step1_case_dedup_report.json`

## Summary Metrics
- ZIP files processed: 88
- ZIP files with DEMO: 88
- ZIP files with DELETE: 17
- Unique deleted IDs: 85,250
- Raw DEMO rows: 23,992,742
- Rows after delete filtering: 23,934,455
- Rows after primaryid-level dedup: 23,929,056
- Rows after caseid-level dedup: 19,961,057
- Exported latest case rows: 19,961,057

## Validation Checks
- `case_selected` total rows = 19,961,057
- distinct `caseid` = 19,961,057
- distinct `primaryid` = 19,961,057
- duplicate `caseid` groups = 0
- duplicate `primaryid` groups = 0

## Notes
- Legacy AERS quarters required row-length correction due delimiter inconsistencies; corrected rows are tracked per ZIP in `step1_ingest_per_zip.csv`.
- Current output is the authoritative key table for downstream joins on DRUG/REAC tables.
