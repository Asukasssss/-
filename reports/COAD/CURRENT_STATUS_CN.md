# COAD 账号 B 当前进展

本轮已完成全部35基因的细胞来源扩展：Lee与Pelka各34个基因可评估，GSTT2均缺独立条目，按不可评估保留；共85,820行。HDC/GSTA4与旧版逐字段一致，新增32个可测基因；ADSS2按官方旧名ADSS匹配。全35展示与完整分层：`results/COAD/06_EXTERNAL/20260920T081200Z_cell_expression35_v2/README_CN.md`。无新P/q，974/39及9/21/5不变，外部代谢关系验证仍未完成。


最新细胞来源批次已完成：Lee63,689细胞/23患者、Pelka371,223细胞/62患者，按组织、富集、技术及患者分开汇总HDC/GSTA4，共4,904行。HDC的肿瘤肥大细胞定位由Pelka支持；Lee肿瘤该类仅3细胞/2患者，不能作为第二队列稳健验证。GSTA4在两研究的基质与上皮中较明显；没有验证CAMP代谢轴或估计组织贡献比例。结果：`results/COAD/06_EXTERNAL/20260920T072000Z_cell_expression_v1/README_CN.md`。实际PRMT7两探针均匹配两个参考转录本，论文V1/V2对应仍未解决；见05_FUNCTION/20260920T073000Z_prmt7_probe_v1。方向表补了AQP9原研究条件，未新增独立研究。974/39及9/21/5不变；下一步优先HDC肥大细胞定位、GSTA4上皮/基质区分，空间验证与纯度仍未完成。

最新协变量批次：M1年龄+性别、M2再加分期各652/652条可计算；q<0.05分别0/0条。均n=33，M0_CC精确复用原统计，未因协变量缺失减样本。M3原部位类别1/2人，按运行前支持规则不可评估。完整974行及原39行字段保持不变。结果：`results/COAD/04_ROBUSTNESS/20260919T142500Z_covariates_v1/README_CN.md`。纯度/技术批次及外部验证仍未完成。

最新功能交付：用户提供的35基因/39关系功能核查包已完成文件核对、本端40条来源标识/摘要复核、6条正文选读与真实高精度回接。原39行和974行每个原字段保持不变；完整池74行获得基因级文献注释，逐关系范围仍为39。研究安排保留9优先、21保留、5暂挂；05/07整阶段均PARTIAL。

补充BCAT2条件性移植瘤证据及正常BCAA饮食单独敲低的阴性结果，不能概括为普适抑瘤；UCKL1功能不等于尿苷介导。完整接受报告及解释：`results/COAD/05_FUNCTION/20260919T134848Z_function_review_accepted_v1/README_CN.md`；当前表为该目录gene_function_review.tsv和*_annotated_exact.tsv。原Excel未改，原作者包保存在imported_delivery作来源快照；实际GitHub交付见发布记录。

当前候选口径（用户在查看结果后指定）：按主关联 P<0.05 选择 **39条探索性优先关系，涉及35个基因、23个代谢特征**；25条正相关、14条负相关。全部标记“名义P显著、未通过FDR”；原q及完整974行保留，不对子集重新校正。两种敏感性均P<0.05且同向的有33条，仅作描述，不再缩小本批39条名单。下一步优先核查35个基因的适用功能证据，其他候选仍保留取证资格。

探索性名单与八节说明：`results/COAD/07_INTEGRATION/20260919T125421Z_nominal_priority_v1/README_CN.md`。这是看过结果后按用户要求改变优先核查标准，未改变下面的首轮统计结论；最终证据整合仍为PARTIAL。

首轮未调整数值批次已完成：33 名经 GEO individual 核对的 I–IV 期患者，预设 674 条关系中 652 条可计算，新的关联 q<0.05 为 0 条（最小 q=0.326）。446 个基因可做33对肿瘤—正常RNA比较，314 个 q<0.05，其中146上升、168下降。RNA差异与肿瘤内部相关不是同一问题，不据此认定全部候选无效。

全部37个作者Tumor（含3腺瘤、1个0期）敏感性及可用性敏感性均为0条 q<0.05。可用性当前版本为 v2：437条相同输入复用主分析，215条样本减少项保留原统计；所有974候选行保留。652条均完成逐一剔除患者检查。M1/M2协变量敏感性已完成，M3不可评估；纯度/技术批次及外部验证尚未完成。

患者结果：`results/COAD/03_PATIENT/20260919T122034Z_patient_v1/README_CN.md`；可用性修正：`results/COAD/04_ROBUSTNESS/20260919T122940Z_availability_reuse_v2/README_CN.md`。正式顺序与状态见 coordination/stages/COAD.tsv。数值已核查；GitHub发布事实见 coordination/publications/COAD.json，分支可查不等于合并main。

本轮从 159 条冻结效应记录出发，对其中原 q<0.05 的 73 条重新建立数据库证据。未复用原基因映射、BRCA 候选或其他癌种分析。冻结统计不变，86 条非显著背景保留。

| 已完成环节 | 本轮结果 |
|---|---|
| 名称与标识初查 | 73 条全部记录去向；63 条原 KEGG 名称相容，另 10 条有身份、暴露背景或新编号核查事项 |
| KEGG 反应候选发现 | 60 个特征，693 条特征—人类基因候选关系；尚不等于直接底物证据 |
| 人类催化反应注释核对 | 482 条候选匹配到 reviewed 人类蛋白的具体 Rhea 反应；其中 58 条存在身份暂挂 |
| 结构化转运反应检索 | 43 个特征，281 条关系；其中 31 条存在身份暂挂 |
| 两类证据去重合并 | 763 条有反应注释支持；其中 674 条无当前特定身份暂挂，涉及 59 个特征、458 个基因 |

**下表为独立生化映射覆盖；患者结果见上方最新批次。当前仍不是最终靶点排名，没有身份暂挂不等于完成原始鉴定。**

# 会改变下一步分析的发现

1. glutamate 的 KEGG 是 L 型，原 HMDB 是 D 型：保留原记录与两个假设，暂不把相关基因放入正式检验家族。
2. ribulose/xylulose 5-phosphate 是合并特征：不能拆成两个独立测量。
3. cysteine-glutathione disulfide 原 KEGG 栏是反应编号 R00900；Mucate 的原 C01807 本次未返回；另为 ophthalmate 查到拟议 C21016。这些提议均未覆盖原注释。
4. Disulfiram 等五条需先核对原始注释或暴露背景，不能因药物有靶点就直接连入代谢关系。
5. 人类反应核对保留了 beta-alanine—UPB1 等有直接反应注释的组合；游离 proline—P4HA1 没有匹配到相同反应，留在待审区，不能把蛋白残基修饰当作游离代谢物转换。
6. 转运检索补充了 taurine—SLC6A6、creatinine 对应转运关系等，避免仅凭酶反应缺失就淘汰特征。

# 全部 73 条的当前覆盖

| 原特征 | 原 g | 原 q | 催化注释基因数 | 转运注释基因数 | 身份状态 |
|---|---:|---:|---:|---:|---|
| beta-alanine | 2.226 | 3.399e-10 | 5 | 8 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| serotonin (5HT) | -2.577 | 5.92e-10 | 2 | 10 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| glucose | -2.468 | 7.024e-10 | 17 | 16 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| Disulfiram | -2.024 | 1.467e-08 | 0 | 0 | HOLD_EXPOSURE_OR_ANNOTATION_REVIEW |
| hypotaurine | 1.358 | 3.777e-08 | 5 | 4 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| proline | 1.832 | 3.777e-08 | 6 | 11 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| taurine | 1.768 | 3.777e-08 | 5 | 9 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| N-acetylaspartate (NAA) | -1.609 | 4.01e-08 | 4 | 2 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| S-adenosylmethionine (SAM) | 1.587 | 9.785e-08 | 33 | 1 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| 3PG | -1.622 | 2.5e-07 | 8 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| cytidine 5'-diphosphocholine | 1.203 | 7.433e-07 | 4 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| creatinine | -1.216 | 2.154e-06 | 0 | 4 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| ribulose/xylulose 5-phosphate | -1.389 | 2.185e-06 | 9 | 0 | HOLD_COMBINED_FEATURE |
| 2-phosphoglycerate | -1.445 | 2.185e-06 | 9 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| 4-methyl-2-oxopentanoate | -1.475 | 2.228e-06 | 4 | 3 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| allantoin | -1.133 | 2.759e-06 | 1 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| urate | -1.382 | 3.037e-06 | 1 | 13 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| acetylcholine | -1.264 | 3.037e-06 | 2 | 5 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| adenosine | -1.434 | 3.684e-06 | 8 | 7 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| histamine | -1.217 | 4.589e-06 | 3 | 5 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| 3-hydroxy-3-methylglutarate | 1.283 | 4.645e-06 | 0 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| citrate | -1.082 | 4.645e-06 | 5 | 7 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| glycylglycine | 1.094 | 5.325e-06 | 0 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| N-acetylglucosamine 6-phosphate | 1.206 | 1.513e-05 | 4 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| Acetohydroxamate | -0.521 | 1.513e-05 | 0 | 0 | HOLD_EXPOSURE_OR_ANNOTATION_REVIEW |
| glutamine | -1.185 | 1.532e-05 | 15 | 19 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| hypoxanthine | 1.134 | 7.858e-05 | 4 | 3 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| guanosine | -1.106 | 7.945e-05 | 4 | 6 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| 5-hydroxylysine | 1.102 | 0.0001065 | 1 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| adenine | 1.006 | 0.0001719 | 3 | 4 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| Mucate | -1.073 | 0.0001719 | 0 | 0 | HOLD_ORIGINAL_ENTRY_NOT_RETURNED |
| kynurenine | 1.036 | 0.0002183 | 3 | 2 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| UDP-N-acetylglucosamine | 0.939 | 0.0002908 | 6 | 7 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| cis-aconitate | -0.835 | 0.0003078 | 1 | 1 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| ophthalmate | -0.854 | 0.0003127 | 1 | 0 | PROVISIONAL_NEW_KEGG_ID |
| isocitrate | -0.940 | 0.0007013 | 5 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| N-acetylneuraminate | -0.914 | 0.0009277 | 2 | 1 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| adenylosuccinate | 0.723 | 0.001036 | 3 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| inosine 5'-monophosphate (IMP) | -0.889 | 0.001036 | 22 | 1 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| 5-methylthioadenosine (MTA) | 0.794 | 0.001107 | 5 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| uridine | -0.929 | 0.001112 | 9 | 6 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| glutathione, reduced (GSH) | 0.724 | 0.00138 | 43 | 4 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| urea | -0.857 | 0.001917 | 3 | 6 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| glutamate | 0.753 | 0.002053 | 49 | 31 | HOLD_ID_CONFLICT |
| triethanolamine | -0.763 | 0.002053 | 0 | 0 | HOLD_EXPOSURE_OR_ANNOTATION_REVIEW |
| dihydroxyacetone phosphate (DHAP) | -0.921 | 0.00209 | 8 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| uracil | 0.735 | 0.003175 | 3 | 4 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| inosine | -0.467 | 0.003175 | 8 | 5 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| guanine | -0.590 | 0.003401 | 5 | 4 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| N-acetylmethionine | 0.760 | 0.00352 | 0 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| creatine | -0.803 | 0.00511 | 6 | 3 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| phosphoenolpyruvate (PEP) | -0.682 | 0.005825 | 9 | 2 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| guanidinoacetate | -0.523 | 0.007479 | 2 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| 1-methylnicotinamide | 0.737 | 0.01043 | 1 | 2 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| 1-methyladenosine | 0.745 | 0.01043 | 0 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| UDP-glucuronate | 0.514 | 0.01043 | 21 | 3 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| Phthalate | -0.608 | 0.01043 | 0 | 0 | HOLD_EXPOSURE_OR_ANNOTATION_REVIEW |
| glycerophosphorylcholine (GPC) | -0.672 | 0.01131 | 8 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| 2-Deoxyglucose 6-phosphate | 0.660 | 0.01268 | 0 | 0 | HOLD_EXPOSURE_OR_ANNOTATION_REVIEW |
| 3-hydroxyproline | 0.612 | 0.01287 | 1 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| asparagine | 0.691 | 0.01635 | 4 | 10 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| N-Acetylglucosamine 1-phosphate | 0.652 | 0.02332 | 2 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| 1-methylhistamine | -0.597 | 0.02587 | 2 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| azelate (nonanedioate) | -0.383 | 0.02771 | 0 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| 6-phosphogluconate | -0.642 | 0.02912 | 4 | 0 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| cysteine-glutathione disulfide | 0.592 | 0.02954 | 0 | 0 | HOLD_REACTION_ID_IN_COMPOUND_FIELD |
| carnosine | -0.552 | 0.03032 | 3 | 5 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| spermidine | 0.306 | 0.03827 | 4 | 8 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| leucine | 0.508 | 0.0399 | 4 | 17 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| tyrosine | 0.544 | 0.0399 | 6 | 7 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| isoleucine | 0.568 | 0.04461 | 4 | 7 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| uridine 5'-diphosphate (UDP) | 0.401 | 0.04585 | 56 | 3 | ORIGINAL_KEGG_NAME_COMPATIBLE |
| cytidine | -0.485 | 0.04585 | 7 | 5 | ORIGINAL_KEGG_NAME_COMPATIBLE |

# 未完成与接续条件

- 作者注释来源、标本对应、GEO患者单位及配对已核对；3个官方旧符号恢复后，预设家族652条可计算。不能把身份暂挂因来源相符就解除。
- 年龄、性别、分期和部位可用，协变量分析尚未运行；纯度与技术批次仍待核实。
- 部分HMDB来源访问受限；底物特异性、复合体角色和不完整数据库覆盖仍待补查。无命中不当阴性。
- 35基因定向功能核查已接入，完整候选功能覆盖、DepMap及外部独立资料仍未完成；新文献不改变原关联q。

# 文件与来源

- 当前特征目录：current_catalog_v0_1/feature_catalog.json；完整候选（含暂挂和待审）：current_catalog_v0_1/candidate_catalog.json。
- 身份与 KEGG 反应明细：fresh_mapping_v0_1/；人类蛋白逐条核对：human_reaction_check_v0_1/；转运反应：transport_check_v0_1/。每条证据含 URL，来源清单含访问时间和哈希。
- 脚本：code/coad/ 下对应脚本。正式七阶段索引：coordination/stages/COAD.tsv；GitHub 发布事实：coordination/publications/COAD.json。早期八阶段 JSON 包为历史格式，不是另一套共同规范。
- 原数据 SHA-256：`1bb1d57480a0fbc2185d11f7598e67e7443aef9e40ea8a036c6da4e5503a4597`。
- 核查通过：原效应/原 q/原标识逐值未改；73 条全覆盖；唯一特征—基因键无重复；UPB1 与 SLC6A6 阳性关系、P4HA1 底物区分和 glutamate 暂挂规则通过。

公开来源：[KEGG](https://www.kegg.jp/kegg/rest/keggapi.html)、[HMDB 谷氨酸条目](https://hmdb.ca/metabolites/HMDB0003339)、[UniProt](https://www.uniprot.org/help/api)、[Rhea](https://www.rhea-db.org/help/download)。
