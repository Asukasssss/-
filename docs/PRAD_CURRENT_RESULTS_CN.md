# PRAD分析入口

## 当前交付：2026-09-22第一轮

- **[完整报告与文件入口](../results/PRAD/07_INTEGRATION/20260922T142000Z_discovery_v1/README_CN.md)**
- [全部101基因易读比较表](../results/PRAD/07_INTEGRATION/20260922T142000Z_discovery_v1/全部101基因_易读比较表.tsv)
- [完整161条直接关系](../results/PRAD/07_INTEGRATION/20260922T142000Z_discovery_v1/candidate_relations_integrated.tsv)
- [单细胞来源与图](../results/PRAD/06_EXTERNAL/20260922T142000Z_discovery_v1/README_CN.md)

43对、361项代谢物：69项P<0.05，9项q<0.05。48项映射161条关系、101基因。肿瘤关联18条名义显著、0条FDR支持；配对RNA21基因FDR支持。单细胞一项24供者肿瘤研究覆盖100基因，91基因可解释来源排名；GPX1缺测。

仍待补：21项特征的直接映射/身份、第二独立肿瘤单细胞研究。无独立患者关系或功能验证，不称机制闭环。服务器原矩阵与其他癌种未修改。

用户于2026-09-22指定本任务负责PRAD，参考本项目BRCA分析。独立分支：`analysis/prad-discovery-20260922`；服务器根：`results/collaborative/PRAD/B/`。

主线遵循[来源SOP](https://github.com/Asukasssss/-/blob/6dd0a73/docs/CAMP_DISCOVERY_TO_CELL_SOURCE_SOP_CN.md)，BRCA参考提交`f18a215dc8a269e4f16c55011633d557f069f631`。本任务不修改BRCA或其他癌种结果，不改变共同标准。

依次登记：输入与作者病例配对；全部保留代谢物配对发现及敏感性；P<0.05工作池直接映射；全关系肿瘤RNA关联与全候选RNA背景；全候选单细胞表达来源；整合。缺测保留状态，不当作阴性。

运行版本及实际完成状态见`coordination/stages/PRAD.tsv`。原始矩阵、完整标本连接、病例测量均仅在server165。
