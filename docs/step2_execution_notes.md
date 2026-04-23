# Step2 Execution Notes

Run window: 2026-04-11 18:21:38 to 2026-04-12 08:25:56
Script: `scripts/step2_build_case_drug_reaction.py`

## Inputs
- Step1 keys: `outputs/step1/faers_latest_case_keys.csv`
- Raw ZIPs: 88 quarters (2004Q1-2025Q4)

## Outputs
- `outputs/step2/case_drug_filtered_std.csv`
- `outputs/step2/case_reaction_filtered_std.csv`
- `outputs/step2/step2_report.json`
- `outputs/step2/step2_quarter_stats.csv`
- `outputs/step2/step2_quarter_stats_clean.csv`
- `outputs/step2/reaction_unresolved_top.csv`
- `outputs/step2/faers_step2.sqlite`

## Core results (from step2_report.json)
- `meddra_loaded`: false
- `latest_keys_rows`: 19,961,057
- `case_drug_rows`: 72,212,791
- `case_reaction_rows`: 59,395,334
- `case_drug_distinct_primaryid`: 19,742,296
- `case_reaction_distinct_primaryid`: 19,961,005

Reaction mapping source counts:
- `NO_MEDDRA_DICTIONARY`: 33,376,683
- `MANUAL_HINT_NO_MEDDRA`: 89

## Important runtime notes
- Pipeline was resumed multiple times; raw `step2_quarter_stats.csv` contains duplicate log fragments for a few quarters.
- Cleaned per-quarter log is provided in `step2_quarter_stats_clean.csv`.
  - raw rows: 182
  - clean rows: 176 (= 88 quarters x 2 tables)

## Validation highlights
- Quarter coverage complete to 2025Q4 for both DRUG and REAC.
- Step1 key consistency (sampled modulus check):
  - sampled case_drug missing key pairs: 0
  - sampled case_reaction missing key pairs: 0

## Pending
- MedDRA dictionary mapping (PT/LLT -> PT code/name) pending local dictionary approval.
