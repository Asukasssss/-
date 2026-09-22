# PDAC 当前进度：作者细胞身份整改 v4

## 本轮问题

接回作者恶性/正常标签，并更正此前把三队列统一为上皮/导管后丢失的身份解释。

## 输入与范围

原687基因、148777细胞、715关系保留。GSE242230采用cell_type_specific；其他两队列已核对现有GEO元数据只有导管大类。原始计数和完整逐细胞映射留服务器。

## 实际结果

- GSE242230：宽类Malignant的11133细胞中，细分为5660 Classical、5123 Basal和350 Normal Epithelial。现10783个作者恶性与350个作者正常上皮分开计算；350/11133=3.14%。
- GSE263733与GSE278688：保留Ductal (unresolved)，不能将肿瘤组织来源直接等同恶性身份。作者论文做过恶性分析，但当前没有接回完整逐细胞判定。
- 三队列严格同首位类别120基因，同首位且各自保持率≥80%为72基因（旧158）；当前250直接基因中为30（旧65）。身份类别不等价导致的变化不称生物学证据消失。
- GSE242230拆分后146基因的首位跨宽谱系改变；完整标签变化见source_identity_changes.tsv。正常上皮350细胞中291细胞分布于6个合格来源标签，恶性10783细胞中10759细胞分布于23个合格来源标签，两组不是匹配对照。
- 配对代谢物仍307项/51项P<0.05/0项q<0.05；当前RNA仍247项可评估/91项P<0.05/40项q<0.05。357直接、358条件关系及所有关联P/q不变。
- ABAT在GSE242230的首位为Myeloid；完整候选来源与比较表已更新，不以ABAT等少数候选决定规则。

## 新手解释

“作者注释恶性”保留作者的证据来源；“身份未定”是尚未得到逐细胞判定，既不称正常也不称恶性。细分正常上皮标签优先于冲突的宽标签。旧29个上皮/导管稳定来源基因不能写成29个三队列恶性来源复现。

## 限制/反证

本轮未自行运行CNV或做机制分析。两队列身份连接仍待解决；临床身份去重与外部代谢关系验证未因本次整改完成。正常上皮350细胞分布及每类≥3来源标签/每标签≥20细胞阈值会限制可评估性；不可评估不当作零。

## 当前决定

GSE242230已完成作者身份接回、重新计算、独立核对和全候选图。06身份总任务仍PARTIAL；另两队列NEEDS_REVIEW。v4来源解释优先，旧数值不覆盖。PR #2仍草稿、未合并。

## 下一步

取得其他两队列可连接到细胞条码的作者恶性/正常判定后，另建版本接入；不自动用表达高低或肿瘤样本标签替代。

## 复现命令

服务器run_author_identity_v4_server.py；本地deliver_author_identity_v4.py。保留完整参数、输入哈希、逐细胞连接核查和旧新比较。

- [身份整改结果](../../results/PDAC/06_EXTERNAL/20260922T131800Z_author_identity_v4/README_CN.md)
- [来源图](../../results/PDAC/06_EXTERNAL/20260922T131800Z_author_identity_v4/figures/README_CN.md)
- [关系比较](../../results/PDAC/07_INTEGRATION/20260922T131900Z_author_identity_v4/candidate_relations_integrated.tsv)
- [基因比较](../../results/PDAC/07_INTEGRATION/20260922T131900Z_author_identity_v4/candidate_genes_integrated.tsv)
- [恶性与正常表达汇总](../../results/PDAC/07_INTEGRATION/20260922T131900Z_author_identity_v4/author_malignant_normal_gene_profiles.tsv)
- [完整旧新来源比较](../../results/PDAC/07_INTEGRATION/20260922T131900Z_author_identity_v4/source_identity_changes.tsv)
