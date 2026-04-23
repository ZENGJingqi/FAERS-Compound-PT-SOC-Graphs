# Graph Releases

## Released graph archives

1. `step6_compound_pt_soc_core_ge10.sqlite.gz`
2. `step6_compound_pt_soc_core_ge20.sqlite.gz`
3. `step6_compound_pt_soc_core_ge30.sqlite.gz`

## Shared schema

### Node types
- `compound`
- `pt`
- `soc`

### Edge types
- `compound_has_pt_ps`
- `compound_has_pt_ss`
- `pt_belongs_to_primary_soc`

## Edge weights

For `compound -> pt` edges, the weight is `n_reports`, defined as the number of FAERS reports where the standardized compound and standardized PT co-occur under the given role (`PS` or `SS`).

## Recommended default release

The recommended default release is `ge20`.

## Counts

| Graph version | compound | pt | soc | PS edges | SS edges |
|---|---:|---:|---:|---:|---:|
| Full cleaned data | 3732 | 21674 | 27 | 1807206 | 2092482 |
| n_reports >= 10 | 2006 | 9917 | 27 | 339013 | 459064 |
| n_reports >= 20 | 1707 | 7274 | 27 | 202917 | 279990 |
| n_reports >= 30 | 1564 | 6035 | 27 | 148522 | 207842 |
