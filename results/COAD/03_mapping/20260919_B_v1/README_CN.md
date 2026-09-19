# 03_mapping 直接生化关系映射

## 基本信息

COAD / B / 03_mapping / 20260919_B_v1

批次状态：DONE（独立生化候选初稿；保留支持、暂挂、待审和不可评估）。

## 输入与分析范围

独立生化候选初稿；保留支持、暂挂、待审和不可评估。

输入：reports/COAD/current_catalog_v0_1/candidate_catalog.json

## 实际结果

记录数：974；单位：特征—人类基因组合。

{"all_candidate_pairs": 974, "reaction_supported_pairs": 763, "supported_without_current_identity_hold": 674, "features_without_current_identity_hold": 59, "genes_without_current_identity_hold": 458, "record_status_counts": {"DONE": 674, "IN_PROGRESS": 288, "NOT_EVALUABLE": 12}}

## 结果含义

DONE 表示本批指定的描述/数据库初查完成；不等于原始鉴定、样本关联或功能验证。

## 校验

公共字段、顺序、记录数、唯一键、输入及输出哈希由 tools/validate_stage_results.py 校验。

## 未完成与限制

- 原冻结效应/原 q 不改；身份与底物范围未完全核实。
- 患者矩阵与凭据不在此包。
- server165 本机别名无法解析，RNA/样本分析未做。

## 下一步

核对作者身份与样本覆盖，补查待审关系后确定下一版分析范围。

## 发布状态

ACCESS_BLOCKED：本地结果已准备，尚无可用写入凭据，不能称为已上传。后续发布状态查 coordination/publications/COAD.json。
