# BRCA单细胞结果集中交付

## 本轮问题
将单细胞数值、来源稳定性、分型、热图、点图与全156基因UMAP集中呈现。
## 输入与范围
7个已发布运行目录，源提交见analysis_spec。未读取患者矩阵。
## 实际结果
完整文件数见validation，逐文件哈希见single_cell_artifact_catalog.tsv。
## 新手解释
先看完整图册，再用页码索引找基因，最后对照来源表的供者覆盖与表达差距。
## 限制/反证
仅整理入口，无新统计或新图；表达来源不代表功能或代谢物来源。
## 当前决定
全156均保留，单细胞停在来源；重点PNG不是唯一候选。
## 下一步
从统一候选表讨论有限问题，不自动新增机制分析。
## 复现
python code/brca_single_cell_catalog_v1.py；再次发布用新RUN_ID。

[单细胞图表总入口](../../../../docs/BRCA_SINGLE_CELL_RESULTS_CN.md)
