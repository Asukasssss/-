# COAD 02_MAPPING

## 本轮问题

将现有 COAD 结果对齐 main 已发布的七阶段规范；不重算或另做候选筛选。

## 输入与范围

原冻结效应表及本账号独立建立的数据库证据；具体来源和哈希见 source_manifest.tsv。结果表只是同一来源的格式视图。

## 实际结果

73 条显著特征身份初查；974 条候选组合，763 条有反应注释支持，其中 674 条当前无特定身份暂挂，涉及 59 个特征和 458 个基因。98 条身份暂挂、190 条反应特异性待查、12 条 reviewed 人类条目不可评估。身份明细见 reports/COAD/fresh_mapping_v0_1/identity_review.json。

## 新手解释

原 g 是标准化效应，不是倍数；原 q 不是新的关联 q。数据库反应注释不等于 COAD 患者支持或功能成立。

## 限制/反证

原始身份、部分底物特异性与复合体角色未全部核实；患者身份和矩阵未读取。server165 别名无法解析。results.tsv 的 NA 及 status/reason 保留缺项，不表示阴性。

## 当前决定

效应描述批次完成；映射整体 PARTIAL，保留全部待审和不可评估组合。统一展示不改变历史数值或 A/BRCA 结果。

## 下一步

获得 server165 正确连接信息后核对作者注释、样本身份和 RNA 覆盖；功能取证可并行。

## 复现命令

校验已生成包：`python code/coad/export_shared_delivery.py --run-id 20260919T112700Z_shared_delivery_v1 --stage 02_MAPPING --validate-only`。重新导出时去掉 --validate-only 并使用新的 UTC run-id；程序拒绝覆盖已有目录。发布状态查 coordination/publications/COAD.json。
