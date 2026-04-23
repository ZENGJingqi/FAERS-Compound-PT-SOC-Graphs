# Graph Releases

This file summarizes the three public FAERS compound-PT-SOC core graph releases.

| Graph version | Minimum `n_reports` threshold | Compound nodes | PT nodes | SOC nodes | PS edges | SS edges |
|---|---:|---:|---:|---:|---:|---:|
| ge10 | 10 | 2006 | 9917 | 27 | 339013 | 459064 |
| ge20 | 20 | 1707 | 7274 | 27 | 202917 | 279990 |
| ge30 | 30 | 1564 | 6035 | 27 | 148522 | 207842 |

Recommended default: `ge20`.

- `ge10` keeps more coverage but includes more low-frequency edges.
- `ge20` is the default release because it balances coverage and noise control.
- `ge30` is the strictest release and keeps only stronger compound-PT links.
