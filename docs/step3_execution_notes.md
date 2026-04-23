# Step3 Execution Notes

Run date: 2026-04-12
Script: `scripts/step3_map_drug_rxnorm_drugbank.py`

## Objective
Map standardized drug terms from Step2 to:
1. RxNorm (`RXNCONSO.RRF`)
2. DrugBank IDs (via RxNorm `SAB=DRUGBANK` links)
3. DrugBank InChI structures

## Inputs
- `outputs/step2/faers_step2.sqlite` (`case_drug` table)
- `RxNorm_full_03022026/rrf/RXNCONSO.RRF`
- `drugbank_5.1.15/XML_dbid_dname.csv`
- `drugbank_5.1.15/XML_dbid_inchi.csv`

## Outputs
- `outputs/step3/faers_step3.sqlite`
- `outputs/step3/step3_report.json`
- `outputs/step3/drug_term_final.csv`
- `outputs/step3/unmapped_drug_terms_top.csv`
- `outputs/step3/term_stage_counts.csv`

## Final coverage (after enhanced-rule full rerun)
From `outputs/step3/step3_report.json` (latest run at 2026-04-12 15:19):
- Unique drug terms: 702,719
- Terms mapped to RxNorm: 226,517 (32.23%)
- Terms mapped to DrugBank ID: 196,419 (27.95%)
- Terms mapped to InChI: 183,162 (26.06%)
- Terms manually excluded as non-drug: 17,305

Row-weighted (using term frequencies from case_drug):
- Drug rows total: 72,212,745
- Drug rows mapped to RxNorm: 68,723,303 (95.17%)
- Drug rows mapped to DrugBank ID: 64,460,834 (89.27%)
- Drug rows mapped to InChI: 53,858,885 (74.58%)
- Drug rows manually excluded as non-drug: 1,094,195

## Key implementation notes
- Multi-stage matching: `EXACT -> CLEANED -> PAREN -> SPLIT`.
- Added separator normalization for backslash-combined ingredient strings (`A\B`), which significantly improved row-level coverage.
- Mapping is deterministic and reproducible from local files only.
- Added incremental script: `scripts/step3_unmapped_incremental_curation.py`.
- Added `RXNREL` bridge matching (`RXNREL.RRF`) to infer DrugBank-linked RXCUIs via ingredient/tradename/form relations, not only direct `RXNCONSO` links.
- Added enhanced normalization rules: salt suffix stripping, component-wise cleanup for slash combinations, and stricter unresolved fallback.
- Incremental full pass over all unresolved terms (`599,813` candidates) added:
  - `26,638` new InChI-mapped terms
  - `17,146` non-drug exclusions
- Manual curation round on current high-frequency unresolved top300:
  - File: `resources/step3_manual/corrections_batch2_top300_manual.csv`
  - Decisions: `MAP_TO_TERM=290`, `EXCLUDE_NON_DRUG=10`
  - Applied result: `map_success=107`, `exclude=10`, `unresolved=183`
- Batch progress logs persisted under:
  - `resources/step3_manual/batches/`
  - `resources/step3_manual/decisions_top5000.csv`
  - `resources/step3_manual/corrections_master.csv`

## Latest enhanced run summary
- Run mode: `incremental_unmapped_all`
- Selected unresolved terms in this run: `555,912`
- Rule-based updates:
  - `map_success=53,659`
  - `exclude=1`
  - `unresolved=502,252`

## Pending
- No SMILES yet (current structure endpoint is InChI).
- MedDRA-independent path is complete; reaction MedDRA coding will be added after dictionary approval.
