# Step5 执行记录（MedDRA 反应标准化）

更新时间：2026-04-13

## 1. 目标
- 在不重跑全量 Step2 的前提下，基于现有 `case_reaction` 做 MedDRA v29.0 标准化。
- 将 `pt_norm` 统一映射到标准 PT 代码与术语，构建可直接用于后续分析的标准化反应底座。

## 2. 输入与环境
- Step2 数据库：`outputs/step2/faers_step2.sqlite`
- MedDRA 字典：`D:\博士文件\TCMMKG\data\AEMS_FDA不良反应数据\MedDRA\MedDRA_29_0_English\MedAscii`
- MedDRA 版本：`29.0$English$$$$`

## 3. 执行过程

### 3.1 增量标准化（脚本）
- 脚本：`scripts/step5_standardize_reaction_meddra.py`
- 逻辑：
  - 提取 `case_reaction` 中去重后的 `pt_norm`（24,675 个）；
  - 映射策略：`PT_EXACT` -> `LLT_TO_PT` -> `MANUAL_CODE_HINT`；
  - 生成表：
    - `reaction_term_meddra_map`
    - `reaction_unresolved_top`
    - `reaction_base_meddra`（仅 mapped）
    - `pt_case_counts_meddra`

### 3.2 人工补丁（小规模未映射）
- 原始未映射：19 个词条（2,391 条记录）。
- 新增人工覆盖文件：`resources/step5_manual/reaction_term_manual_overrides_2026-04-13.csv`
- 脚本：`scripts/step5_apply_manual_reaction_overrides.py`
- 结果：
  - 覆盖候选 18，实际应用 18；
  - 额外恢复映射记录：1,679 条；
  - 最终仅剩 1 个未映射词条（`VAGINAL CUFF`，1 条记录）。

### 3.3 重建标准化后的结构-反应关系
- 脚本：`scripts/step5_build_inchikey_pt_meddra.py`
- 输入：
  - `outputs/step4/faers_step4.sqlite`（`case_inchikey_psss`）
  - `outputs/step5/faers_step5.sqlite`（`reaction_base_meddra`）
- 输出：
  - MedDRA 标准化的 `InChIKey-PT(PS/SS)` 关系表及建模矩阵。

## 4. 关键结果

### 4.1 反应标准化结果（最终）
- `reaction_rows_total`: 59,395,334
- `reaction_rows_mapped`: 59,395,333
- `reaction_rows_unmapped`: 1
- `reaction_row_coverage_mapped`: 99.999998%
- `term_total`: 24,675
- `term_mapped`: 24,674
- `term_unmapped`: 1
- `term_coverage_mapped`: 99.995947%
- `pt_distinct_mapped`: 22,643
- `primaryid_distinct_mapped`: 19,961,005

### 4.2 标准化后 InChIKey-PT（PS/SS）
- `inchikey_case_counts_psss_meddra_rows`: 3,758
- `inchikey_pt_counts_psss_meddra_rows`: 2,692,460
- `inchikey_pt_psss_meddra_model_rows`（`n_cases>=20`）: 394,254
- `inchikey_index_meddra_rows`: 1,748
- `pt_index_meddra_rows`: 7,833
- `matrix_edges_rows`: 394,254

## 5. 产物文件
- `outputs/step5/faers_step5.sqlite`
- `outputs/step5/step5_report.json`
- `outputs/step5/step5_manual_override_report.json`
- `outputs/step5/reaction_term_meddra_map.csv`
- `outputs/step5/reaction_unresolved_top.csv`
- `outputs/step5/pt_case_counts_meddra.csv`
- `outputs/step5/step5_structure_report.json`
- `outputs/step5/inchikey_pt_counts_psss_meddra.csv`
- `outputs/step5/inchikey_pt_counts_psss_meddra_min20.csv`
- `outputs/step5/inchikey_index_meddra.csv`
- `outputs/step5/pt_index_meddra.csv`
- `outputs/step5/inchikey_pt_matrix_psss_meddra_min20.csv`

## 6. 当前状态结论
- 反应标准化已达到“可发表级别”的清洁度（仅 1 条未映射记录）。
- 后续可直接进入：
  - 信号统计（ROR/PRR/IC/EBGM）；
  - DDI 结构对风险表构建；
  - 外部验证（JADER/说明书证据）。
