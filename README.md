# FAERS Preprocessing and Core Heterograph Package

This repository-ready package contains only data collection notes, preprocessing code, curated documentation, reference tables, and three pruned FAERS core graphs. It does not include any downstream modeling or prediction results.

## Scope

The package documents how raw FAERS data were transformed into standardized compound-PT-SOC graphs using:

- FAERS / AERS quarterly ASCII archives from 2004Q1 to 2025Q4
- RxNorm `RxNorm_full_03022026`
- DrugBank `5.1.15`
- MedDRA `29.0` English and Chinese

## Main outputs

- `docs/`: process notes and execution notes
- `figures/`: English workflow and graph summary figures
- `tables/`: English summary tables and cleaned reference tables
- `graphs/`: three pruned graph archives (`ge10`, `ge20`, `ge30`) compressed as `.sqlite.gz`
- `code/`: final preprocessing and graph-building scripts

## Graph definition

Each core graph keeps only:

- node types: `compound`, `pt`, `soc`
- edge types: `compound_has_pt_ps`, `compound_has_pt_ss`, `pt_belongs_to_primary_soc`

`compound -> pt` edges are weighted by `n_reports`, defined as the number of FAERS reports where the compound and PT co-occur under the given role (`PS` or `SS`).

## Graph versions

The package provides three pruned graph versions:

- `ge10`: keep `compound -> pt` edges with `n_reports >= 10`
- `ge20`: keep `compound -> pt` edges with `n_reports >= 20`
- `ge30`: keep `compound -> pt` edges with `n_reports >= 30`

The recommended default graph is `ge20`, which gives the best tradeoff between coverage and noise control.

## Files required for reproduction

- `tables/data_collection_versions.csv`
- `tables/graph_summary.csv`
- `tables/stage_changes.csv`
- `tables/data_overview.csv`
- `tables/all_compounds_basic_info.xlsx`
- `tables/all_reactions_basic_info_bilingual.xlsx`
- `tables/all_soc_basic_info_bilingual.xlsx`
- `tables/all_compound_pt_ps_associations_bilingual.csv.gz`
- `tables/all_compound_pt_ss_associations_bilingual.csv.gz`
- `code/`
- `graphs/`

## Notes about graph archives

Graph archives are compressed as `.sqlite.gz` to make GitHub publication feasible. Use standard gzip tools to decompress them:

```bash
gzip -d compound_pt_role_core_ge20.sqlite.gz
```
