# 各癌种阶段结果与发布规范 v1.0

用户要求每阶段上传 GitHub，并保持跨癌种顺序、格式、字段一致。机器定义以 `config/stage_result_schema_v1.json` 为准。

## 固定阶段顺序

01 效应记录 → 02 身份核对 → 03 直接生化关系 → 04 样本/患者关联 → 05 功能/反证 → 06 外部证据 → 07 候选分层 → 08 重点深入。

这是展示与归档顺序，不要求串行。功能和外部取证可与患者分析并行。未做阶段在索引中保持同一位置并记 NOT_RUN。

## 固定目录与文件

入口：`results/<CANCER>/stage_index.json`。每批：`results/<CANCER>/<stage_id>/<batch_id>/`，新版本另建目录，不覆盖已发布批次。

每批固定四文件：README_CN.md（schema 中八标题顺序）、summary.json（固定字段顺序）、records.json（固定公共前缀、result 保存详细结果）、manifest.json（前三文件的字节数和 SHA256，不包含自身）。

COAD 的 reports/COAD 保留分析明细和来源；标准入口供跨癌种读取。BRCA 历史 reference 快照保持只读，由 A 在后续交付中按规范生成；B 不擅自改写 A 的结果。

## 字段、顺序和含义

- records 的公共字段严格遵循 record_prefix。缺标识用 null，缺数值不能填 0。
- original_effect、original_q 保留冻结表精确字符串。新 RNA/相关 P、q 放在 result 中独立字段，附 analysis_type 和 test_family。
- 标准记录按 dataset、feature_name、metabolite_key、gene_symbol、human_gene_id 升序；Unicode 字符序、null 视为空字符串。报告可以另外展示按 q 排序的榜单。
- 不因同名合并特征。关系键为队列×实际特征×基因×版本；多个反应/来源不构成重复检验。
- summary 中说明实际单位，不能混淆代谢特征、关系和基因数量。
- record_status 仅说明本批指定工作的完成度；DONE 不等于身份鉴定、患者支持或功能成立。细分原因放 reason_code。
- 使用 DONE、IN_PROGRESS、NOT_RUN、NOT_EVALUABLE、ACCESS_BLOCKED。阶段进度与批次分开：初查批次 DONE 不表示身份核对或全部映射已完成。

## 每阶段上传

1. 生成四文件，更新本癌种索引、交接说明。
2. 运行 `python tools/validate_stage_results.py --cancer <CANCER>`、`python tools/check_repository.py`；核对新内容可共享、哈希、数量、字段、排序和唯一键。
3. 每阶段/可交付批次建立独立提交，推送工作分支，通过 PR 交接。不能等项目全部做完再上传，不强制推送，不自动合并共享规则改动。
4. 记录 branch、远程 commit、PR URL 与核实时间。PUSHED、PR_OPEN、MERGED 分开；不能把本地完成称为上传完成。
5. 认证或网络失败记 ACCESS_BLOCKED，并向用户说明缺少的条件。

批次 publication 是生成时快照；之后的发布事实追加到 `coordination/publications/<CANCER>.json`，引用提交与批次，不为更新发布状态重写已发布结果。最新状态需核对远程 refs/PR。

只上传代码、配置、说明、数据库来源与允许共享的小型汇总；不上传 runtime 缓存、患者级数值、源矩阵、完整 sample mapping 或凭据。
