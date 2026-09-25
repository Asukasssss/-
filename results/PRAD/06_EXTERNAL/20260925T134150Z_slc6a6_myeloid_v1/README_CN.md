# SLC6A6 髓系亚群直接细胞比较

本轮问题：作者是否有可追溯的髓系亚群注释；同亚群肿瘤与癌旁细胞的SLC6A6表达如何？用户后续明确选择直接细胞比较，故本轮不进行供者配对或等权。

输入与范围：原PRAD24固定CELLxGENE版本，髓系共2743细胞，7个minor类别和13个subset。原有字段与作者元数据字典、20221004注释脚本和20230714更新代码一致，包含mature_DC更名mDC_regs、CLEC10A/FCGR3A/CLEC9A/IRF4等更新。当前沿用作者v2注释，不重跑聚类，不以单一标志基因重新定类。来源代码固定提交 f48298f64d5f4eb02709b21969e95a0014cc8c09，来源/哈希见source_manifest.tsv。

实际结果：
- Macrophages: cancer n=583,adjacent n=436; means=0.172388/0.206842; detection=19.4%/25.2%; P=0.0610131
- Mast_cells: cancer n=64,adjacent n=51; means=0.105625/0.0943846; detection=6.2%/7.8%; P=0.783751
- Monocytes: cancer n=525,adjacent n=340; means=0.347785/0.363875; detection=35.4%/34.1%; P=0.939326
- cDC1: cancer n=24,adjacent n=35; means=0.348972/0.262959; detection=41.7%/40.0%; P=0.561187
- cDC2: cancer n=348,adjacent n=285; means=0.223524/0.222397; detection=27.6%/29.1%; P=0.81841
- mDC_regs: cancer n=13,adjacent n=8; means=0.690257/0.621718; detection=61.5%/50.0%; P=0.850448
- pDC: cancer n=14,adjacent n=17; means=0/0.13478; detection=0.0%/23.5%; P=0.0610137

新手解释：每类细胞分别合并全部肿瘤细胞与全部癌旁细胞，展示平均表达和检测比例，Mann-Whitney检验比较分布。P<0.05按细胞层面标记，q值保留。不同患者贡献细胞数不等，P不是患者层面的证据。minor主比较与subset探索比较分别列族，不只挑显著项。

限制/反证：样本仍来自同一研究；癌旁不是健康人。作者标签可追溯，但不是独立金标准验证，特别是少量DC亚群。标志基因图仅作表达一致性检查，不是新注释或功能证据。作者泛化ontology中的部分DC仅标为myeloid cell；本轮使用更细的作者标签而非将ontology强行细化。此前每侧20细胞、至少8供者配对标准下，细分亚群覆盖不足；本轮分析单位依用户要求改为细胞，不声称克服了供者覆盖限制，也不覆盖既往供者分析。表达不等于转运活性。

当前决定：仅用于直接表达对照与假设生成，保留所有结果及细胞数。原有供者配对报告不修改。

下一步：若进行患者层面的验证，仍需更多配对供者/独立数据。

复现：python prad_slc6a6_myeloid_cells_v1.py --out NEW_RUN --source-run ORIGINAL_RUN --code-commit COMMIT --author-commit AUTHOR_COMMIT。源矩阵、细胞值、细胞身份留server165；公开汇总与图。
