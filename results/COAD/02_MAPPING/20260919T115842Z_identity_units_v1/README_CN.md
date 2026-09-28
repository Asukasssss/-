# COAD 患者单位与基因符号核对

## 本轮问题

核对 37 个作者 Tumor 标本的患者单位、病变构成和 22 个 RNA 未匹配基因 ID。核对批次 DONE；映射与关联阶段整体尚未完成。

## 输入与范围

使用 GEO GSE89076 的样本元数据、作者 MasterMapping、既有 RNA 行名及 HGNC 官方记录。以 RNAID 精确连接 GEO accession，并核对 GEO title 与 MetabID；患者身份来自明确的 individual 字段，不从数字或排序猜配对。完整临床连接表及 GEO 样本元数据只在服务器运行目录保存。

## 实际结果

- GEO 的 80 个样本来自 39 个明确 individual；其中两人各有两个不同阶段的肿瘤取样。作者当前映射排除了这两人的肿瘤取样，纳入的 37 个 Tumor 标本对应 37 个不同 individual，并有 37 个按 individual 精确对应的正常对照。
- 37 个 Tumor 中：腺瘤 3、0 期 1、I 期 8、II 期 8、III 期 9、IV 期 8。年龄、性别、部位、分期字段完整。
- 3 个官方旧符号唯一恢复：GBA1→RNA 行 GBA，SLC35D4→TMEM241，SLC60A2→MFSD4B。HGNC entrez_id 一致，旧符号反查均只命中同一 HGNC 记录。
- hsa:102724197 已通过 NCBI 确认为 LOC102724197，描述为 inactive glutathione hydrolase 2；本矩阵没有其精确行，不借用 GGT 家族其他行。
- 其余未恢复项目保留，具名缺项不直接称为平台未测或基因不表达。完整 974 条中 938 条现在有 RNA 对应。当前有反应注释支持且无特定身份暂挂的 674 条中，652 条有 RNA，涉及 446/458 个基因 ID。

## 新手解释

独立单位核对在作者提交的 GEO 注释层面成立，尚无基因型指纹验证。作者 Tumor 标签包含腺瘤和 0 期，因此不能把全部 37 个直接称为同质的浸润癌组。官方改名恢复的是数据对应，不是新的生化证据或功能证明。

## 限制/反证

化合物身份暂挂均未解除；反应特异性和复合体角色仍待核对。旧符号唯一性不证明微阵列探针无交叉杂交。缺少匹配的 readthrough 不拆分成组成基因。674 条只是本版本的暂定关联家族，不代表最终直接机制关系。

## 当前决定

在看关联结果前固定主分析为 I–IV 期 33 个肿瘤标本；全部 37 个作者 Tumor 作为病变构成敏感性分析。674 条关系均列入预设家族并记录 22 条 RNA 缺项；BH 对本家族实际可计算项执行，计划数与实际数分开。原 CAMP 效应/q 不改，新的分析另标版本。

## 下一步

运行肿瘤内部 Spearman 关联、配对 RNA 对照和预设敏感性分析。缺项、身份暂挂以及功能证据继续整理，相关显著不是功能候选唯一门槛。

## 复现命令

脚本：`code/coad/fetch_gene_identity_review.py`、`code/coad/resolve_units_and_symbols.py`。服务器脚本必须在新建的 COAD/B 独占目录执行，参数 `--project-root <根目录> --run-dir <新目录> --code-commit <脚本提交>`。来源及 SHA256 见 source_manifest.tsv、gene_source_registry.json；公开汇总哈希见 validation.json。不得同步 private_clinical_join.tsv 或 geo_samples.soft。

来源：[GEO GSE89076](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89076)、[原研究](https://pubmed.ncbi.nlm.nih.gov/28847964/)、[HGNC API](https://www.genenames.org/help/rest/)、[NCBI Gene 102724197](https://www.ncbi.nlm.nih.gov/gene/102724197)。
