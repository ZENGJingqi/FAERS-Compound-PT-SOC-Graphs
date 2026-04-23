# FAERS Compound-PT-SOC Graphs

This repository contains the Step1-Step6 preprocessing workflow, three final FAERS core graph releases, and three lightweight reference tables used to organize standardized compound-PT-SOC relationships. It does not include any downstream modeling or prediction results.

## Workflow Overview

![Workflow overview](figures/workflow_overview.png)

## Scope

The released workflow covers only:

1. case deduplication
2. case-drug and case-reaction table construction
3. drug term normalization and DrugBank linkage
4. structure normalization to InChIKey / SMILES
5. reaction normalization to MedDRA PT / primary SOC
6. core graph construction and pruning

## Data Sources and Versions

- FAERS / AERS quarterly ASCII archives: `2004Q1-2025Q4`
- RxNorm: `RxNorm_full_03022026`
- DrugBank: `5.1.15`
- MedDRA English: `29.0`
- MedDRA Chinese: `29.0`

## Repository Structure

- `code/`: final Step1-Step6 preprocessing and graph-building scripts
- `docs/`: English process notes, execution notes, and graph release notes
- `figures/`: workflow and graph summary figures
- `graphs/`: three pruned graph archives compressed as `.sqlite.gz`
- `tables/`: three curated reference tables for compounds, PT terms, and SOC terms

## Expected Local Layout For Reproduction

The scripts assume the repository is placed beside two user-supplied directories:

```text
your_working_directory/
|-- FAERS-Compound-PT-SOC-Graphs/
|-- raw_data/
|   `-- faers_quarterly_archives/
`-- reference_data/
    |-- RxNorm_full_03022026/
    |-- drugbank_5.1.15/
    `-- MedDRA/
```

In other words:

- raw FAERS/AERS ZIP files should be placed in `../raw_data/faers_quarterly_archives/`
- reference datasets should be placed in `../reference_data/`

## Core Graph Schema

Each released core graph keeps only:

### Node types

- `compound`
- `pt`
- `soc`

### Edge types

- `compound_has_pt_ps`
- `compound_has_pt_ss`
- `pt_belongs_to_primary_soc`

`compound -> pt` edges are weighted by `n_reports`, defined as the number of FAERS reports where the standardized compound and standardized PT co-occur under the given role (`PS` or `SS`).

## Graph Versions

- `step6_compound_pt_soc_core_ge10.sqlite.gz`: keep `compound -> pt` edges with `n_reports >= 10`
- `step6_compound_pt_soc_core_ge20.sqlite.gz`: keep `compound -> pt` edges with `n_reports >= 20`
- `step6_compound_pt_soc_core_ge30.sqlite.gz`: keep `compound -> pt` edges with `n_reports >= 30`

The recommended default graph is `ge20`, because it provides the best balance between coverage and noise control.

## Graph Summary

![Graph comparison](figures/graph_comparison.png)

## Reference Tables Included In This Public Repository

The repository now includes three lightweight reference tables:

- `tables/all_compounds_basic_info.xlsx`
- `tables/all_reactions_basic_info_bilingual.xlsx`
- `tables/all_soc_basic_info_bilingual.xlsx`

These files provide the cleaned compound identifiers, bilingual PT names, and bilingual SOC names used by the released graphs.

## Included Files

- `README.md`
- `code/`
- `docs/`
- `figures/`
- `graphs/`
- `tables/`
- `release_manifest.json`

See [docs/GRAPH_RELEASES.md](docs/GRAPH_RELEASES.md) for the released graph summary.

## Files Intentionally Not Uploaded

This public repository still does **not** include the full PS/SS association matrices:

- full PS compound-PT association table
- full SS compound-PT association table

## Decompression

```bash
gzip -d step6_compound_pt_soc_core_ge20.sqlite.gz
```
