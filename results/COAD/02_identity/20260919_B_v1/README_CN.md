# 02_identity 代谢物身份核对

## 基本信息

COAD / B / 02_identity / 20260919_B_v1

批次状态：DONE（数据库名称与标识初查；不代表源鉴定全部完成）。

## 输入与分析范围

数据库名称与标识初查；不代表源鉴定全部完成。

输入：reports/COAD/fresh_mapping_v0_1/identity_review.json

## 实际结果

记录数：73；单位：显著代谢特征。

{"record_status_counts": {"DONE": 63, "IN_PROGRESS": 10}}

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
