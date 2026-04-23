# Reaction Standardization Strategy (FAERS)

## Why this is critical
Reaction terms are the core labels for downstream ADR analysis. Noise in PT strings directly degrades signal quality and model validity.

## What we borrowed from prior work
1. `WangLabCSU/faers`:
- Standardize REAC/INDI terms against MedDRA hierarchy.
- Use curated typo-to-code corrections before dictionary matching.
- Prefer PT exact match, then LLT to PT roll-up.

2. `faersdbstats`:
- Map reaction PT to standard MedDRA concepts (`vocabulary_id='MedDRA'`).
- Build standardized case-outcome and drug-outcome tables for signal statistics.

## Implemented in this project (Step2)
### A. Deterministic text normalization
For each REAC `pt`:
- Unicode normalize (NFKC)
- trim
- collapse spaces
- uppercase
- normalize punctuation/hyphen variants

Output: `pt_norm`

### B. Curated manual correction hints
- Integrated curated reaction term -> MedDRA code hints from `faers` package `clean_reac_pt()` logic.
- Output field: `pt_manual_code_hint`

### C. Optional MedDRA dictionary mapping (deferred for now)
If MedDRA dictionary exists locally:
1. PT exact: `pt_norm` -> PT code
2. LLT fallback: `pt_norm` -> LLT -> PT code
3. Manual code hint fallback: hint code -> PT name

Outputs:
- `pt_meddra_code`
- `pt_meddra_pt`
- `pt_mapping_source`

## Current status without MedDRA dictionary
- Dictionary-dependent mapping was intentionally skipped (`meddra_loaded=false`).
- Reaction standardization still yields stable `pt_norm` labels.
- This is safe for interim counting/co-occurrence work.
- Final MedDRA concept-level harmonization will be executed after dictionary approval.

## Data-quality implications
- Interim phase: lexical standardization + dedup + case filtering already remove major noise.
- Final publication phase: add MedDRA PT/LLT mapping and optionally SOC/HLGT roll-up for clinical hierarchy analysis.
