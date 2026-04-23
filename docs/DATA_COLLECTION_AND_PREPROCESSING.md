# Data Collection, Preprocessing, and Core Graph Construction

## 1. Data collection

The preprocessing workflow used the following local datasets:

1. FAERS / AERS quarterly ASCII archives from 2004Q1 to 2025Q4
2. RxNorm `RxNorm_full_03022026`
3. DrugBank `5.1.15`
4. MedDRA `29.0` English
5. MedDRA `29.0` Chinese

## 1.1 Expected reproduction layout

To run the public scripts as documented, place the cloned repository beside `raw_data/` and `reference_data/`:

```text
your_working_directory/
├─ FAERS-Compound-PT-SOC-Graphs/
├─ raw_data/
│  └─ faers_quarterly_archives/
└─ reference_data/
   ├─ RxNorm_full_03022026/
   ├─ drugbank_5.1.15/
   └─ MedDRA/
```

This is the directory convention assumed by the default script arguments in `code/`.

## 2. Step1-Step6 workflow

### Step1. Case deduplication
- remove duplicate and superseded FAERS case versions
- retain the latest valid case record
- result: 19,961,057 latest valid cases

### Step2. Case-drug and case-reaction tables
- build normalized `case_drug` and `case_reaction` tables
- result: 72,212,791 case_drug rows and 59,395,334 case_reaction rows

### Step3. Drug term normalization and DrugBank linkage
- normalize raw drug terms
- map to RxNorm / DrugBank-compatible names
- manually curate high-frequency unmapped terms
- anchor structure linkage using DrugBank InChI
- result: 702,719 unique drug terms, 183,162 terms mapped to InChI, 53,858,885 drug rows mapped to InChI

### Step4. Structure normalization
- standardize linked structures
- retain `InChIKey`, `inchi_std`, and `smiles_std`
- use `InChIKey` as the unique compound identifier
- result: 4,523 unique InChIKey, with 3,758 observed in PS/SS reports

### Step5. Reaction normalization
- normalize reaction terms to MedDRA Preferred Terms (PT)
- assign MedDRA primary System Organ Class (SOC)
- add bilingual PT and SOC names using MedDRA English and Chinese
- result: 59,395,333 reaction rows mapped to MedDRA PT, with 22,643 distinct PT terms

### Step6. Core graph construction
The final core graph keeps only:
- nodes: `compound`, `pt`, `soc`
- edges: `compound_has_pt_ps`, `compound_has_pt_ss`, `pt_belongs_to_primary_soc`
- edge weight: `n_reports` for `compound -> pt`

## 3. Full cleaned reference data retained locally

The local working project retains:

- all compounds: 3,732 rows
- all PT nodes: 21,674 rows
- all SOC nodes: 27 rows
- all PS compound-PT associations: 1,807,206 rows
- all SS compound-PT associations: 2,092,482 rows

## 4. Public graph versions

The public repository includes three pruned graph releases:

- `ge10`
- `ge20`
- `ge30`

These versions are defined by the minimum `n_reports` threshold applied to `compound -> pt` edges.

## 5. Recommended graph

The recommended default graph is `ge20`, because it offers the best balance between lower noise and broader coverage.
