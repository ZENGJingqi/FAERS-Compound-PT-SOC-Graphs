# Step4 执行记录（InChIKey-Reaction 数据集）

更新时间：2026-04-12

## 1. 目标
- 在 Step3 已有 `InChI` 映射基础上，构建可用于建模的一阶段结构-反应数据：
  - 以 `InChIKey` 作为结构唯一标识；
  - 统一输出标准 `SMILES` / `InChI_std`；
  - 将 `InChIKey` 与 Step2 的 `case-reaction(pt_norm)` 建立可统计关系；
  - 暂不做 MedDRA 字典标准化（待本地字典到位后进入下一步）。

## 2. 执行逻辑
- 输入：
  - `outputs/step2/faers_step2.sqlite`
  - `outputs/step3/faers_step3.sqlite`
- 核心处理：
  - 仅保留 Step3 中有 `InChI` 且非 `MANUAL_EXCLUDE` 的药物词条；
  - 用 RDKit 进行标准化（cleanup / fragment parent / uncharge / canonical tautomer）；
  - 生成 `InChIKey`、`SMILES_std`、`InChI_std`；
  - 按 `PS/SS` 角色构建 `InChIKey-PT` 计数；
  - 生成建模矩阵索引（`inchikey_id`, `pt_id`, `n_cases`）。

## 3. 本轮实际执行说明
- 先清空了 `outputs/step4` 旧内容后重建。
- 为保证可恢复性，新增了两个脚本：
  - `scripts/step4_resume_post_inchi.py`（从中间表断点续跑）
  - `scripts/step4_finalize_pt_export.py`（仅收尾 PT 聚合+导出）
- 最终以 `step4_finalize_pt_export.py` 完成收尾并输出完整文件。

## 4. 关键结果（来自 step4_report.json）
- `term_inchi_raw_rows`: 183,162
- `raw_inchi_norm_ok`: 4,788 / 4,800（失败 12）
- `inchikey_meta_rows`（唯一 InChIKey）: 4,523
- `case_inchikey_any_rows`: 43,415,079
- `case_inchikey_psss_rows`: 19,600,684
- `inchikey_pt_psss_rows`: 2,788,221
- `inchikey_pt_psss_model_rows`（min_pair_cases=20）: 398,908
- `inchikey_index_rows`: 1,745
- `pt_index_rows`: 8,254
- `matrix_edges_rows`: 398,908
- `include_any`: false（本轮只保留 PS/SS 主分析链路）

## 5. 产物目录
- `outputs/step4/faers_step4.sqlite`
- `outputs/step4/step4_report.json`
- `outputs/step4/inchi_normalization_map.csv`
- `outputs/step4/term_inchikey_map.csv`
- `outputs/step4/inchikey_meta.csv`
- `outputs/step4/inchikey_case_counts_any.csv`
- `outputs/step4/inchikey_case_counts_psss.csv`
- `outputs/step4/inchikey_pt_counts_psss_min20.csv`
- `outputs/step4/inchikey_pt_counts_psss_top20000.csv`
- `outputs/step4/inchikey_pt_counts_any_top20000.csv`
- `outputs/step4/inchikey_index.csv`
- `outputs/step4/pt_index.csv`
- `outputs/step4/inchikey_pt_matrix_psss_min20.csv`

## 6. 下一步
- 等 MedDRA 本地字典到位后，对 Step2/Step4 使用的 `pt_norm` 做 MedDRA 标准化（LLT/PT层级）；
- 在 MedDRA 标准化版本上重建 `InChIKey-PT`，得到第一阶段可发表的数据底稿。
