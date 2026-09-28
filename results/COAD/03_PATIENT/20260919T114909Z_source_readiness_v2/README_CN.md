# COAD 作者数据与候选覆盖核对

## 本轮问题

服务器连接恢复后，核对 COAD 原始作者文件、标本对应和候选数据覆盖，并检查冻结注释是否忠实来自作者文件。本批来源核对 DONE；03_PATIENT 整体 PARTIAL，患者关联 NOT_RUN。

## 输入与范围

在 server165 的独占目录 `results/collaborative/COAD/B/20260919T114909Z_source_readiness_v2/` 执行。输入包括作者 MasterMapping、PreprocessedData_COAD.xlsx、GSE89076 的处理后 RNA 矩阵，以及本账号已有的冻结效应、身份初查和候选目录。精确路径及 SHA256 见 source_manifest.tsv；脚本和版本见 analysis_spec.json。

仅对作者明确映射的 37 个肿瘤标本汇总候选两端的有限数值数量，不下载患者级矩阵，不推断配对，不新增转换、填补、相关检验或 BH 校正。

## 实际结果

| 核查项 | 结果 |
|---|---|
| 作者 COAD 映射 | 76 行：37 肿瘤、39 正常；三个标本标识字段均无缺失、无重复 |
| 肿瘤映射覆盖 | 37 个标本全部对应到代谢和 RNA 矩阵 |
| 作者处理后代谢矩阵 | 肿瘤 160×37，正常 160×39；159 条冻结特征全部找到 |
| RNA 矩阵 | 21,634 个基因行、80 个样本列；行列标识唯一 |
| 原注释来源 | 73 条显著特征均唯一匹配作者 H_name，KEGG/HMDB 与冻结表逐值一致 |
| 候选基因 | 644 个不同人类基因 ID；622 个按精确符号找到，21 个具名符号未匹配，另 1 个 ID 尚无解析符号 |
| 全部候选组合 | 974 条：935 条两端均有 37 个有限数值；36 条 RNA 符号未匹配；3 条基因符号未解析 |
| 674 条当前无特定身份暂挂且有反应支持的组合 | 649 条两端均有数据，25 条 RNA 符号未匹配；仍不是正式检验家族 |
| 作者 sampleanno | 550 行（275 TUMOR、275 NORMAL），仅 SAMPLE_NAME/GROUP 两列，无显式患者 ID 字段 |

完整结果见 results.tsv；分项见 feature_coverage.tsv、gene_coverage.tsv、identity_source_review.tsv。前三列数值效应与 q 仅在 original_effect/original_q 中保留；本批 effect/p_value/q_value 均为 NA。

## 新手解释

935 条表示目前“有两端测量数据”，不是 935 条显著关联。37 是作者映射的肿瘤标本数，不能称为已确认的 37 名独立患者。原始 sampleanno 的 550 行也不能直接当作本轮有两种组学的分析样本数。

注释与作者原表一致，说明先前的 glutamate L/D 型冲突、合并峰和反应编号问题并非本次导入产生；一致性不能证明化学身份正确，也不能解除暂挂。

## 限制/反证

患者独立性、肿瘤—正常配对、协变量、RNA 全部预处理步骤及作者 data 页的原始检出含义仍未确认。不按编号推断患者或配对。当前精确 RNA 符号未匹配不等于基因不表达；别名、改名和平台覆盖需另外核查。974 条全部保留，ready_for_patient_testing 均为 false。

服务器 v1 核对保留于其原目录；本 v2 将一个符号未解析的基因 ID（影响 3 条组合）从“RNA 符号未匹配”中单列。总覆盖数与冻结统计未变，v1/v2 不构成独立证据。

## 当前决定

连接阻塞已解除。保留原身份暂挂和候选范围，不因缺失符号删除条目；不提前锁定正式检验家族。02_MAPPING 通过阶段索引引用本批身份来源核对，不复制为另一批独立结果。

## 下一步

核对原研究的患者单位、样本选取及可用协变量；为 21 个未匹配符号和 1 个未解析 ID 核查权威标识，记录唯一映射或无法评估原因。随后固定肿瘤内部分析范围、版本和检验家族。COAD 功能取证可并行。

## 复现命令

本地验证汇总包：`python tools/validate_coad_source_readiness.py --run-id 20260919T114909Z_source_readiness_v2`。

服务器复跑必须新建独占目录，复制脚本与同版公开输入，建立内容为新 run-id 的 `.running` 文件，再执行 `python3 <新目录>/audit_server_sources.py --project-root <服务器项目根> --run-dir <新目录> --code-commit <脚本提交号>`。脚本拒绝已完成目录；不复写本批。

来源读取前后 SHA256 一致、标本连接、记录数、唯一键、冻结效应/q 字符串及导出哈希已校验，见 validation.json。发布事实见 coordination/publications/COAD.json；阶段上传与 main 合并分开记录。
