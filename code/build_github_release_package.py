from __future__ import annotations

import gzip
import hashlib
import json
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DOCS = ROOT / 'docs'
SOURCE_SCRIPTS = ROOT / 'scripts'
SOURCE_SUMMARIES = ROOT / 'outputs' / 'summaries'
SOURCE_REFERENCE = ROOT / 'outputs' / 'reference_tables'
SOURCE_GRAPHS = ROOT / 'outputs' / 'step6_core_graphs'
TARGET = ROOT / 'release' / 'FAERS-Compound-PT-SOC-Graphs'

plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


def sha256sum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def ensure_clean_dir(path: Path) -> None:
    if not path.exists():
        path.mkdir(parents=True, exist_ok=True)
        return
    for child in path.iterdir():
        if child.name == '.git':
            continue
        if child.is_dir():
            shutil.rmtree(child)
        else:
            child.unlink()


def gzip_file(src: Path, dst: Path) -> None:
    with src.open('rb') as fin, gzip.open(dst, 'wb', compresslevel=6) as fout:
        shutil.copyfileobj(fin, fout, length=1024 * 1024)


def plot_english_workflow(path: Path) -> None:
    fig, ax = plt.subplots(figsize=(14, 5.6))
    ax.axis('off')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    labels = [
        'Step1\nCase deduplication',
        'Step2\nCase-drug and\ncase-reaction tables',
        'Step3\nDrug name normalization\nRxNorm/DrugBank',
        'Step4\nStructure normalization\nInChIKey/SMILES',
        'Step5\nReaction normalization\nMedDRA PT/SOC',
        'Step6\nPruned core graphs\nge10/ge20/ge30',
    ]
    colors = ['#B23A48', '#464F5F', '#6F4C9B', '#2C5AA0', '#0D7587', '#6B8E23']
    positions = [
        (0.16, 0.78),
        (0.50, 0.78),
        (0.84, 0.78),
        (0.16, 0.24),
        (0.50, 0.24),
        (0.84, 0.24),
    ]
    widths = [0.20, 0.24, 0.24, 0.22, 0.22, 0.24]
    height = 0.21
    for (x, y), label, color in zip(positions, labels, colors):
        width = widths[len(ax.patches)]
        rect = plt.matplotlib.patches.FancyBboxPatch(
            (x - width / 2, y - height / 2),
            width,
            height,
            boxstyle='round,pad=0.02,rounding_size=0.03',
            fc=color,
            ec='black',
            lw=1.8,
            transform=ax.transAxes,
        )
        ax.add_patch(rect)
        ax.text(
            x,
            y,
            label,
            ha='center',
            va='center',
            color='white',
            fontsize=18,
            fontweight='bold',
            transform=ax.transAxes,
            linespacing=1.15,
        )

    arrow_kw = dict(arrowstyle='->', lw=2.2, color='black', mutation_scale=24)
    ax.annotate('', xy=(0.36, 0.78), xytext=(0.27, 0.78), xycoords=ax.transAxes, arrowprops=arrow_kw)
    ax.annotate('', xy=(0.70, 0.78), xytext=(0.62, 0.78), xycoords=ax.transAxes, arrowprops=arrow_kw)
    ax.annotate('', xy=(0.24, 0.34), xytext=(0.76, 0.68), xycoords=ax.transAxes, arrowprops=arrow_kw)
    ax.annotate('', xy=(0.39, 0.24), xytext=(0.28, 0.24), xycoords=ax.transAxes, arrowprops=arrow_kw)
    ax.annotate('', xy=(0.72, 0.24), xytext=(0.61, 0.24), xycoords=ax.transAxes, arrowprops=arrow_kw)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
    fig.savefig(path, dpi=300, bbox_inches='tight', pad_inches=0.03)
    plt.close(fig)


def plot_graph_comparison(path: Path, comp_df: pd.DataFrame, ref_summary: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    plot_df = pd.concat(
        [
            pd.DataFrame([
                {
                    'label': 'Full cleaned data',
                    'compound_nodes': ref_summary['compound_rows'],
                    'pt_nodes': ref_summary['reaction_rows'],
                    'soc_nodes': ref_summary['soc_rows'],
                    'ps_edges': ref_summary['ps_association_rows'],
                    'ss_edges': ref_summary['ss_association_rows'],
                }
            ]),
            comp_df.assign(label=comp_df['threshold'].map(lambda x: f'>={x}'))[
                ['label', 'compound_nodes', 'pt_nodes', 'soc_nodes', 'ps_edges', 'ss_edges']
            ],
        ],
        ignore_index=True,
    )
    x = range(len(plot_df))
    labels = plot_df['label'].tolist()

    axes[0].plot(x, plot_df['compound_nodes'], marker='o', color='#B23A48', lw=2, label='compound nodes')
    axes[0].plot(x, plot_df['pt_nodes'], marker='o', color='#2C5AA0', lw=2, label='PT nodes')
    axes[0].plot(x, plot_df['soc_nodes'], marker='o', color='#0D7587', lw=2, label='SOC nodes')
    axes[0].set_xticks(list(x))
    axes[0].set_xticklabels(labels, fontsize=10)
    axes[0].set_ylabel('Node count', fontsize=12, color='black')
    axes[0].tick_params(axis='both', colors='black', labelsize=10)
    axes[0].legend(frameon=False, fontsize=9)
    axes[0].grid(axis='y', linestyle='--', alpha=0.3)

    width = 0.24
    axes[1].bar([i - width / 2 for i in x], plot_df['ps_edges'], width=width, color='#464F5F', label='PS edges')
    axes[1].bar([i + width / 2 for i in x], plot_df['ss_edges'], width=width, color='#6F4C9B', label='SS edges')
    axes[1].set_xticks(list(x))
    axes[1].set_xticklabels(labels, fontsize=10)
    axes[1].set_ylabel('Edge count', fontsize=12, color='black')
    axes[1].tick_params(axis='both', colors='black', labelsize=10)
    axes[1].legend(frameon=False, fontsize=9)
    axes[1].grid(axis='y', linestyle='--', alpha=0.3)

    for ax in axes:
        for s in ax.spines.values():
            s.set_color('black')
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches='tight')
    plt.close(fig)


def main() -> None:
    ensure_clean_dir(TARGET)
    (TARGET / 'docs').mkdir()
    (TARGET / 'figures').mkdir()
    (TARGET / 'graphs').mkdir()
    (TARGET / 'code').mkdir()
    (TARGET / 'tables').mkdir()

    comp_df = pd.read_csv(SOURCE_GRAPHS / 'graph_comparison.csv')
    ref_summary = json.loads((SOURCE_REFERENCE / 'reference_tables_summary.json').read_text(encoding='utf-8'))

    workflow_png = TARGET / 'figures' / 'workflow_overview.png'
    graph_png = TARGET / 'figures' / 'graph_comparison.png'
    plot_english_workflow(workflow_png)
    plot_graph_comparison(graph_png, comp_df, ref_summary)

    for p in SOURCE_SCRIPTS.iterdir():
        if p.is_file() and (p.name.startswith('step') or p.name == 'build_github_release_package.py'):
            shutil.copy2(p, TARGET / 'code' / p.name)

    public_docs = {
        'reaction_standardization_strategy.md',
        'step1_execution_notes.md',
        'step2_execution_notes.md',
        'step3_execution_notes.md',
        'step4_execution_notes.md',
        'step5_execution_notes.md',
        'step6_core_graph_construction.md',
    }
    for p in SOURCE_DOCS.iterdir():
        if p.is_file() and p.name in public_docs:
            shutil.copy2(p, TARGET / 'docs' / p.name)

    for name in [
        'all_compounds_basic_info.xlsx',
        'all_reactions_basic_info_bilingual.xlsx',
        'all_soc_basic_info_bilingual.xlsx',
    ]:
        shutil.copy2(SOURCE_REFERENCE / name, TARGET / 'tables' / name)

    graph_records = []
    for tag in ['ge10', 'ge20', 'ge30']:
        src = SOURCE_GRAPHS / f'core_graph_{tag}' / f'step6_compound_pt_soc_core_{tag}.sqlite'
        dst = TARGET / 'graphs' / f'step6_compound_pt_soc_core_{tag}.sqlite.gz'
        gzip_file(src, dst)
        graph_records.append({
            'graph_version': tag,
            'original_sqlite_bytes': src.stat().st_size,
            'compressed_gzip_bytes': dst.stat().st_size,
            'sha256_gzip': sha256sum(dst),
        })

    readme = """# FAERS Compound-PT-SOC Graphs

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
"""
    (TARGET / 'README.md').write_text(readme, encoding='utf-8')

    methods = """# Data Collection, Preprocessing, and Core Graph Construction

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
|-- FAERS-Compound-PT-SOC-Graphs/
|-- raw_data/
|   `-- faers_quarterly_archives/
`-- reference_data/
    |-- RxNorm_full_03022026/
    |-- drugbank_5.1.15/
    `-- MedDRA/
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

## 6. Public reference tables

The public repository also includes three small reference tables:

1. `tables/all_compounds_basic_info.xlsx`
2. `tables/all_reactions_basic_info_bilingual.xlsx`
3. `tables/all_soc_basic_info_bilingual.xlsx`

These are intended for quick inspection of the released entities without exposing the full large-scale PS/SS association tables.
"""
    (TARGET / 'docs' / 'DATA_COLLECTION_AND_PREPROCESSING.md').write_text(methods, encoding='utf-8')

    graph_release_lines = [
        '# Graph Releases',
        '',
        'This file summarizes the three public FAERS compound-PT-SOC core graph releases.',
        '',
        '| Graph version | Minimum `n_reports` threshold | Compound nodes | PT nodes | SOC nodes | PS edges | SS edges |',
        '|---|---:|---:|---:|---:|---:|---:|',
    ]
    for _, row in comp_df.iterrows():
        graph_release_lines.append(
            f"| ge{int(row['threshold'])} | {int(row['threshold'])} | {int(row['compound_nodes'])} | "
            f"{int(row['pt_nodes'])} | {int(row['soc_nodes'])} | {int(row['ps_edges'])} | {int(row['ss_edges'])} |"
        )
    graph_release_lines.extend([
        '',
        'Recommended default: `ge20`.',
        '',
        '- `ge10` keeps more coverage but includes more low-frequency edges.',
        '- `ge20` is the default release because it balances coverage and noise control.',
        '- `ge30` is the strictest release and keeps only stronger compound-PT links.',
    ])
    (TARGET / 'docs' / 'GRAPH_RELEASES.md').write_text('\n'.join(graph_release_lines) + '\n', encoding='utf-8')

    repo_meta = {
        'repository_name': 'FAERS-Compound-PT-SOC-Graphs',
        'scope': 'Step1-Step6 preprocessing, three final pruned core graphs, and three lightweight reference tables',
        'reference_tables': [
            'all_compounds_basic_info.xlsx',
            'all_reactions_basic_info_bilingual.xlsx',
            'all_soc_basic_info_bilingual.xlsx',
        ],
        'graphs': [
            {
                'graph_version': item['graph_version'],
                'file': f"step6_compound_pt_soc_core_{item['graph_version']}.sqlite.gz",
                'compressed_gzip_bytes': item['compressed_gzip_bytes'],
            }
            for item in graph_records
        ],
    }
    (TARGET / 'release_manifest.json').write_text(json.dumps(repo_meta, ensure_ascii=False, indent=2), encoding='utf-8')
    print(TARGET)


if __name__ == '__main__':
    main()
