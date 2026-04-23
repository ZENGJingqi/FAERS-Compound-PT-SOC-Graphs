# Step6 Core Graph Construction

## Purpose

Step6 constructs the final FAERS core heterographs after compound and reaction standardization. This step keeps only the final data structure required for downstream use and excludes report nodes, compound-compound edges, archetype nodes, and all downstream modeling outputs.

## Final graph schema

### Node types

1. `compound`
2. `pt`
3. `soc`

### Edge types

1. `compound_has_pt_ps`
2. `compound_has_pt_ss`
3. `pt_belongs_to_primary_soc`

### Edge weights

- `compound -> pt` edges are weighted by `n_reports`.
- `n_reports` is the number of FAERS reports in which the standardized compound and standardized PT co-occur under the specified role (`PS` or `SS`).
- `pt -> soc` is a deterministic MedDRA hierarchy mapping.

## Input dependencies

Step6 uses the outputs of:

- Step4: standardized compound structures (`InChIKey`, `smiles_std`, `inchi_std`)
- Step5: standardized MedDRA PT and primary SOC mapping

## Implementation files

1. `step6a_materialize_report_compound_roles_stream.py`
2. `step6b_aggregate_compound_pt_roles.py`
3. `step6c_fix_pt_soc_nodes.py`
4. `step6d_build_pruned_core_graphs.py`
5. `step6e_finalize_bilingual_and_refs.py`

## Final retained graph outputs

### Full standardized graph

Path: `outputs/step6_core_graphs/core_graph_full/step6_compound_pt_soc_core_full.sqlite`

- `compound`: 3,732
- `pt`: 21,674
- `soc`: 27
- `compound_has_pt_ps`: 1,807,206
- `compound_has_pt_ss`: 2,092,482
- `pt_belongs_to_primary_soc`: 21,674

### Pruned graph versions

#### `ge10`
Path: `outputs/step6_core_graphs/core_graph_ge10/step6_compound_pt_soc_core_ge10.sqlite`

- `compound`: 2,006
- `pt`: 9,917
- `soc`: 27
- `compound_has_pt_ps`: 339,013
- `compound_has_pt_ss`: 459,064
- `pt_belongs_to_primary_soc`: 9,917

#### `ge20`
Path: `outputs/step6_core_graphs/core_graph_ge20/step6_compound_pt_soc_core_ge20.sqlite`

- `compound`: 1,707
- `pt`: 7,274
- `soc`: 27
- `compound_has_pt_ps`: 202,917
- `compound_has_pt_ss`: 279,990
- `pt_belongs_to_primary_soc`: 7,274

#### `ge30`
Path: `outputs/step6_core_graphs/core_graph_ge30/step6_compound_pt_soc_core_ge30.sqlite`

- `compound`: 1,564
- `pt`: 6,035
- `soc`: 27
- `compound_has_pt_ps`: 148,522
- `compound_has_pt_ss`: 207,842
- `pt_belongs_to_primary_soc`: 6,035

## Recommendation

The default recommended graph is `ge20` because it provides the best balance between coverage and noise control.
