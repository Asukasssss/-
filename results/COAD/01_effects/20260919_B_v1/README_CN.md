# 01_effects 效应记录整理

## 基本信息

COAD / B / 01_effects / 20260919_B_v1

批次状态：DONE（冻结效应记录描述；不重算效应和 q）。

## 输入与分析范围

冻结效应记录描述；不重算效应和 q。

输入：reference/camp/cancer_effects.tsv

## 实际结果

记录数：159；单位：效应记录。

{"all_effects": 159, "original_q_lt_005": 73, "significant_positive": 33, "significant_negative": 40, "record_status_counts": {"DONE": 159}}

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
