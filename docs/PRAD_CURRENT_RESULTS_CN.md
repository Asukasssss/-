# PRAD分析入口

用户于2026-09-22指定本任务负责PRAD，参考本项目BRCA分析。独立分支：`analysis/prad-discovery-20260922`；服务器根：`results/collaborative/PRAD/B/`。

主线遵循[来源SOP](https://github.com/Asukasssss/-/blob/6dd0a73/docs/CAMP_DISCOVERY_TO_CELL_SOURCE_SOP_CN.md)，BRCA参考提交`f18a215dc8a269e4f16c55011633d557f069f631`。本任务不修改BRCA或其他癌种结果，不改变共同标准。

依次登记：输入与作者病例配对；全部保留代谢物配对发现及敏感性；P<0.05工作池直接映射；全关系肿瘤RNA关联与全候选RNA背景；全候选单细胞表达来源；整合。缺测保留状态，不当作阴性。

运行版本及实际完成状态见`coordination/stages/PRAD.tsv`。原始矩阵、完整标本连接、病例测量均仅在server165。
