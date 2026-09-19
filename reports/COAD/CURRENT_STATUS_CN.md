# COAD 账号 B 当前进展

本轮从 159 条冻结效应记录出发，对其中原 q<0.05 的 73 条重新建立数据库证据。未复用原基因映射、BRCA 候选或其他癌种分析。冻结统计不变，86 条非显著背景保留。

| 已完成环节 | 本轮结果 |
|---|---|
| 名称与标识初查 | 73 条全部记录去向；63 条原 KEGG 名称相容，另 10 条有身份、暴露背景或新编号核查事项 |
| KEGG 反应候选发现 | 60 个特征，693 条特征—人类基因候选关系；尚不等于直接底物证据 |
| 人类催化反应注释核对 | 482 条候选匹配到 reviewed 人类蛋白的具体 Rhea 反应；其中 58 条存在身份暂挂 |
| 结构化转运反应检索 | 43 个特征，281 条关系；其中 31 条存在身份暂挂 |
| 两类证据去重合并 | 763 条有反应注释支持；其中 674 条无当前特定身份暂挂，涉及 59 个特征、458 个基因 |

**这是一版有来源的生化候选目录，不是已完成患者关联的结果或最终靶点排名。没有身份暂挂也不等于完成原始鉴定。所有关系的患者分析均为 NOT_RUN。**

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

- 作者源注释、样本对应、实际 RNA 覆盖与患者独立单位尚未核查；本机 server165 SSH 别名无法解析，等待连接配置。
- 部分 HMDB 页面访问被拒，未完成全面交叉核验；不以网页不返回作为化合物不存在的证据。
- 此映射依赖当前结构化数据库覆盖。无命中不是阴性；KO 未覆盖的人类酶、其他区室表述的转运、复杂底物范围和复合体角色仍需补查。
- 注释中的实验来源不等于我们已阅读论文，也不等于 COAD 功能证明。尚未执行 RNA、患者关联、依赖性或外部验证。
- 接下来先核对源注释和数据覆盖，再锁定版本化关系与检验范围。候选功能取证与患者分析可并行，不把相关显著当唯一门槛。

# 文件与来源

- 当前特征目录：current_catalog_v0_1/feature_catalog.json；完整候选（含暂挂和待审）：current_catalog_v0_1/candidate_catalog.json。
- 身份与 KEGG 反应明细：fresh_mapping_v0_1/；人类蛋白逐条核对：human_reaction_check_v0_1/；转运反应：transport_check_v0_1/。每条证据含 URL，来源清单含访问时间和哈希。
- 脚本：code/coad/ 下对应脚本。新内容目前仅本地，未提交或推送。
- 原数据 SHA-256：`1bb1d57480a0fbc2185d11f7598e67e7443aef9e40ea8a036c6da4e5503a4597`。
- 核查通过：原效应/原 q/原标识逐值未改；73 条全覆盖；唯一特征—基因键无重复；UPB1 与 SLC6A6 阳性关系、P4HA1 底物区分和 glutamate 暂挂规则通过。

公开来源：[KEGG](https://www.kegg.jp/kegg/rest/keggapi.html)、[HMDB 谷氨酸条目](https://hmdb.ca/metabolites/HMDB0003339)、[UniProt](https://www.uniprot.org/help/api)、[Rhea](https://www.rhea-db.org/help/download)。
