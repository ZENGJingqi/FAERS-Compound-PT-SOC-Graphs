# Data Collection, Preprocessing, and Graph Construction

## 1. Data collection

The preprocessing workflow used the following local datasets:

1. FAERS / AERS quarterly ASCII archives from 2004Q1 to 2025Q4
2. RxNorm `RxNorm_full_03022026`
3. DrugBank `5.1.15`
4. MedDRA `29.0` English
5. MedDRA `29.0` Chinese

The source and version information are documented directly in this file for public release.

## 2. Preprocessing workflow

### Step1. Case deduplication

- remove duplicate and superseded FAERS case versions
- retain the latest valid case record

Result:
- latest deduplicated cases: 19,961,057

### Step2. Case-drug and case-reaction tables

- build normalized `case_drug` and `case_reaction` tables after case-level deduplication

Result:
- case_drug rows: 72,212,791
- case_reaction rows: 59,395,334

### Step3. Drug name normalization and DrugBank linkage

- normalize raw drug terms
- map to RxNorm / DrugBank-compatible names
- manually curate high-frequency unmapped terms
- anchor structure linkage using DrugBank InChI

Result:
- unique drug terms: 702,719
- terms mapped to InChI: 183,162
- drug rows mapped to InChI: 53,858,885

### Step4. Structure normalization

- standardize linked structures
- retain `InChIKey`, `inchi_std`, and `smiles_std`
- use `InChIKey` as the unique compound identifier

Result:
- unique InChIKey: 4,523
- InChIKey observed in PS/SS reports: 3,758

### Step5. Reaction normalization

- normalize reaction terms to MedDRA Preferred Terms (PT)
- assign MedDRA primary System Organ Class (SOC)
- add bilingual PT and SOC names using MedDRA English and Chinese

Result:
- reaction rows mapped to MedDRA PT: 59,395,333
- distinct MedDRA PT: 22,643

## 3. Full cleaned reference data retained locally

The local expert package retains the following full cleaned reference data:

- all compounds: 3,732 rows
- all PT nodes: 21,674 rows
- all SOC nodes: 27 rows
- all PS compound-PT associations: 1,807,206 rows
- all SS compound-PT associations: 2,092,482 rows

These full reference tables are not included in the public GitHub repository.

## 4. Core graph construction

The final core graph retains only:

- nodes: `compound`, `pt`, `soc`
- edges: `compound_has_pt_ps`, `compound_has_pt_ss`, `pt_belongs_to_primary_soc`

No report nodes, compound-compound edges, archetype nodes, or downstream model outputs are included.

## 5. Pruned graph versions

Graph versions are defined by the minimum `n_reports` threshold applied to `compound -> pt` edges:

- `ge10`
- `ge20`
- `ge30`

The public repository includes the three graph archives directly under `graphs/`.

## 6. Recommended graph

The recommended default graph is `ge20`, because it offers the best balance between:

- lower noise than `ge10`
- broader coverage than `ge30`
- manageable graph size for downstream use
