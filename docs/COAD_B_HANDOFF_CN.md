# COAD / B 交接

## 当前状态（2026-09-19，首轮患者数值完成）

- DONE：33名明确患者的652条关联（计划674，22条RNA缺项），0条新q<0.05，最小q=0.326；33对RNA有314/446基因q<0.05（146正、168负）。全部974候选行保留，原CAMP统计不变。
- DONE：37人敏感性、652条逐一剔除、可用性子集；两种敏感性均0条q<0.05。可用性v2精确复用437条相同输入主统计和215条样本减少项，再按同一652项重算BH；v1只留追溯。
- 输入/输出哈希、唯一键、公共字段/顺序、独立BH复算及同输入复用已验证。正式路径：`results/COAD/03_PATIENT/20260919T122034Z_patient_v1`、`results/COAD/04_ROBUSTNESS/20260919T122940Z_availability_reuse_v2`。
- PARTIAL：稳健性阶段；协变量调整NOT_RUN，纯度/批次未核实。功能与外部阶段NOT_RUN。无显著关联不当全候选阴性，不改变预设家族来追求显著。
- 下一步：预声明全家族协变量敏感性，批量功能取证与反证核查；患者显著不是唯一门槛。公开文件仅代码/参数/聚合结果，逐患者表留服务器。远程提交与PR事实见 coordination/publications/COAD.json。

## 历史状态（2026-09-19，患者单位与标识确认）

- `20260919T115842Z_identity_units_v1` 核对 DONE：作者纳入的 37 个 Tumor 标本对应 37 个明确 GEO individual，并与正常标本精确配对；原研究两组重复取样未纳入当前 Tumor 集。
- 其中有 3 腺瘤、1 个 0 期，故预先固定 I–IV 期 33 个标本为主分析，全部 37 个为敏感性对照，不更改冻结 CAMP 统计。
- 3 个 HGNC 官方旧符号恢复：GBA1→GBA、SLC35D4→TMEM241、SLC60A2→MFSD4B。完整974关系中938条有RNA；预设674关系中652条有RNA，对应446/458基因。
- 方案 `config/coad_patient_v1.json` 已在统计运行前提交和上传。置换 Spearman、配对 RNA、全部37/作者可用性敏感性分别计算BH；每条主关系提供逐一剔除检查。当前数值计算中，未以此宣称关联完成。

## 历史状态（2026-09-19，服务器来源与覆盖核对）

- 已连接 server165；独占运行 `results/collaborative/COAD/B/20260919T114909Z_source_readiness_v2/` 完成，只同步汇总。
- 作者映射为 37 肿瘤、39 正常标本，映射标识唯一、无缺失，全部对应到作者代谢/RNA矩阵；患者独立性和配对未确认。
- 159 条冻结特征全部找到；73 条显著特征原 KEGG/HMDB 与作者 metanno 一致。身份冲突仍保留，不能因来源一致就解除暂挂。
- 974 条候选中，935 条两端均有 37 个有限值；36 条精确 RNA 符号未匹配，3 条基因符号未解析。涉及 644 个基因 ID，622 个精确符号找到，21 个具名符号未匹配，1 个 ID 未解析。674 条当前支持且无特定身份暂挂组合中 649 条有两端数据。
- 本批来源核对 DONE；03_PATIENT 阶段 PARTIAL，关联 NOT_RUN，正式检验家族尚未锁定。下一步核对患者单位、协变量与基因标识；功能资料可并行。
- 新结果入口：`results/COAD/03_PATIENT/20260919T114909Z_source_readiness_v2/README_CN.md`。v1 的符号缺项分类已在 v2 更正，数值覆盖总量与冻结统计未改；v1 留服务器追溯。

## 历史状态（2026-09-19，连接恢复前的独立生化映射初稿）

- 以 `reports/COAD/CURRENT_STATUS_CN.md` 为最新进展，下面的首批记录保留为历史。
- 从 73 条显著效应记录独立检索 KEGG/HMDB/UniProt/Rhea，未读取旧基因映射或 BRCA 候选。
- DONE：73 条名称/标识初查；693 条 KEGG 反应候选的人类注释核对；2054 个 reviewed 人类 transport 关键词蛋白的结构化转运反应扫描。
- 合并：974 条待评估组合，763 条有催化或转运反应注释支持；其中 674 条当前无特定身份暂挂，涉及 59 个特征、458 个基因。仍属生化候选，非 COAD 关联显著或最终靶点。
- 98 条候选组合因特征身份暂挂（其中 89 条有反应支持），190 条需进一步底物/反应核查，12 条没有关联到 reviewed 人类蛋白条目。各类互斥；无命中不当阴性。
- 影响判断：glutamate 的 KEGG L 型与 HMDB D 型冲突；ribulose/xylulose 是合并特征；R00900 被填入化合物标识栏；Mucate 原 C01807 本次未返回；ophthalmate 新提出 C21016。冻结键和统计未改。
- 五个特征需核查原注释或暴露背景：Disulfiram、Acetohydroxamate、triethanolamine、Phthalate、2-Deoxyglucose 6-phosphate。
- 来源限制：UniProt 实验注释的论文未逐篇复核；复合体角色未逐基因完成；数据库覆盖并非所有可能关系。转运仅纳入结构化 in/out 同一底物且精确 ChEBI/官方 pH 映射的条目。
- ACCESS_BLOCKED：`ssh -o BatchMode=yes -o ConnectTimeout=10 server165 pwd` 返回主机别名无法解析；已向用户请求本机 SSH 连接信息。部分 HMDB 页面/XML 返回 HTTP403。
- NOT_RUN：服务器作者源注释、样本身份、RNA覆盖、患者关联、功能/依赖与外部验证。正式患者检验家族未锁定，ready_for_patient_testing 全为 false。
- 目录：`reports/COAD/current_catalog_v0_1/` 是当前合并目录；其他版本目录保留各阶段明细。公开原始数据库响应缓存位于被忽略的 `runtime/`，无患者数据。每阶段包含脚本、来源 URL/访问时间/哈希。
- 必要核验：冻结哈希及关键字段逐值不变；73 条全覆盖；组合唯一；beta-alanine—UPB1 与 taurine—SLC6A6 阳性对应、游离 proline—P4HA1 未匹配及 glutamate 暂挂断言通过。
- 下一项：连接 server165，核对作者注释和 COAD 样本/RNA覆盖，随后确定有版本和检验范围的分析输入。身份补查与功能资料可并行。
- 已按阶段提交；正式交付入口为 `coordination/stages/COAD.tsv`，与 BRCA 遵循 main 的七阶段规范。历史八阶段 JSON 包仅用于追溯，未修改 A/BRCA 结果、共享映射或冻结表。远程提交与 PR 状态见 `coordination/publications/COAD.json`。

## 首批效应记录梳理（历史）

- 分支：analysis/coad-initial；基线提交：422b75b2f223bb786958f95d6e00d62da98db154。新增内容当前仅本地，尚未提交或推送。
- 用户已指定 B 负责 COAD，并要求从效应记录开始，不直接复用已有映射和候选。
- 本批为本地公开汇总表描述性分析，无服务器运行、无患者级输入。结果：reports/COAD/effect_records_v1/RESULTS_CN.md；参数、哈希、代码版本：同目录 summary.json。
- 输入：reference/camp/cancer_effects.tsv；SHA256：1bb1d57480a0fbc2185d11f7598e67e7443aef9e40ea8a036c6da4e5503a4597。
- 统计对象为冻结效应表中的实际特征记录。本轮没有新假设检验或 BH 校正，未检查独立患者身份。
- DONE：159 条 COAD 记录描述，73 条原 q<0.05，其中 33 正、40 负；86 条不显著背景保留。所有记录均标记 single_dataset。显著记录中 6 条最大原始缺失率>50%，未剔除。
- 判断变化：明确从全部效应记录起步，保留缺失率限制和身份核对需求；尚不形成基因或靶点排名。Disulfiram 名称条目需核查来源，未据此修改冻结注释。
- NOT_RUN：从头核对化合物与直接基因关系、RNA覆盖/作者映射核查、患者关联、功能及外部验证、server165连接核查。
- NOT_EVALUABLE：仅凭本效应表不能确定实际样本量、患者独立性或配对关系；这些信息未在本表提供，不表示服务器一定缺失。
- ACCESS_BLOCKED：本批未遇到。
- 限制：原 HIGH_ID 不是新鉴定；g 不是倍数；通路计数不是富集或通量；原 q 保持不变。
- 下一项工作：B 从 73 条显著特征核对原身份和直接反应来源，独立构建 COAD 映射；不读取旧候选来决定取舍。
- 影响范围：COAD 本地脚本、报告和授权分工登记；未修改 BRCA、公共映射或冻结数据。
