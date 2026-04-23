# Step14 最小主图（compound-pt role-specific core）

日期：2026-04-22

## 目标

按最新确认的口径，构建一版最小、稳定、可直接建模的主异构图：

1. `compound` 只连 `pt`
2. `compound-pt` 必须区分 `PS` 和 `SS`
3. `pt` 再连 `soc`
4. 暂不引入 `report`、`compound-compound`、`compound-archetype`

## 节点

共 `3` 类主节点：

1. `compound`：3732
2. `pt`：10511
3. `soc`：27

说明：

1. `compound` 是在 `PS/SS` 两类角色里实际映射成功并参与主图的成分
2. `pt` 是当前 MedDRA 标准化后、在主图中实际保留的 PT
3. `soc` 是这些 PT 对应的主系统器官分类

## 边

共 `3` 类核心边：

1. `compound_has_pt_ps`：1807206
2. `compound_has_pt_ss`：2092482
3. `pt_belongs_to_primary_soc`：10511

## 边权重

`compound -> pt` 两张边表都保留权重 `n_reports`：

1. `PS` 总权重：36054725
2. `SS` 总权重：50572861

这里的 `n_reports` 表示：

某个成分以 `PS` 或 `SS` 角色出现时，与某个 PT 在同一报告中共同出现的次数。

## 角色覆盖

### PS

1. 源药物记录：19740505
2. 成功映射到结构的源记录：11978551
3. 唯一 `report-compound`：11978551
4. 覆盖报告数：11978550
5. 覆盖成分数：2216

### SS

1. 源药物记录：20966263
2. 成功映射到结构的源记录：14368796
3. 唯一 `report-compound`：9162034
4. 覆盖报告数：4998808
5. 覆盖成分数：3655

## 方法要点

### 1. report-compound 先流式构建

由于直接 SQL join 非常慢，最终采用：

1. 从 `step2.case_drug` 流式读取 `PS`/`SS` 药物记录
2. 用 `step4.term_inchikey_map` 的 `term_norm -> inchikey` 字典在 Python 中映射
3. 写入唯一 `report-compound` 表

这一步稳定跑通了。

### 2. 再按批聚合 compound-pt

将 `report-compound` 按批次与 `step5.reaction_base_meddra` 连接，聚合得到：

1. `edge_compound_pt_ps`
2. `edge_compound_pt_ss`

### 3. pt-soc 直接来自 MedDRA 官方层级

`PT-SOC` 不再依赖中间库，而是直接从：

`MedDRA/MedDRA_29_0_English/MedAscii/mdhier.asc`

中提取 `PT -> primary SOC` 映射，并只保留当前 `node_pt` 中实际出现的 PT。

## 输出

目录：

`D:\博士文件\TCMMKG\data\AEMS_FDA不良反应数据\smiles_adr_project\outputs\step14_heterograph\compound_pt_role_core_v3`

文件：

1. `step14_compound_pt_role_core_v3.sqlite`
2. `report_compound_role_summary.json`
3. `compound_pt_role_summary.json`
4. `node_summary.csv`
5. `edge_summary.csv`
6. `role_summary.csv`

脚本：

1. `step14f2_materialize_report_compound_roles_stream.py`
2. `step14g_aggregate_compound_pt_roles.py`

## 当前判断

这一版主图已经符合当前建模需求：

1. 结构简单
2. 语义清楚
3. 边权重保留了报告频次
4. `PS` 和 `SS` 已明确分开

下一步可以直接围绕这版图设计异构图模型，不需要再回去扩充复杂节点层。
