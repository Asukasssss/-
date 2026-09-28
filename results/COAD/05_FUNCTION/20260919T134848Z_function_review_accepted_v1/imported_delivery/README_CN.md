# COAD 35基因功能核查与研究安排 v1

## 本轮问题
对现有39条名义P<0.05关系涉及的35个基因进行结肠癌/结直肠癌功能核查：
区分真实基因干预、药物或运输证据、其他癌种背景、融合/同名错配、相反结果。
不重新计算患者相关，不把功能证据转换为新的显著q或治疗建议。

## 输入与范围
- 固定仓库：`Asukasssss/-`，COAD分支快照`06921616d5914b6b45c61bdb6c1f259cf1608f84`。
- 原始清单：`results/COAD/07_INTEGRATION/20260919T125421Z_nominal_priority_v1/nominal_candidates.tsv`。
- 原全目录974行、原计划674关系、实际652检验不改。39条原关联全部未通过BH。
- 本轮35基因/39关系的定向检索覆盖完成；并非全部458基因或全部974组合的功能核查。
- 现环境读取的是GitHub汇总及公共论文，未读取server165患者矩阵，也未执行真实974表回接。

## 实际结果
|内容|数量/状态|
|---|---|
|基因检索覆盖|35/35|
|既有名义关系对照|39/39|
|原始来源/排除记录|40条，40个唯一URL|
|有详细来源记录的基因|27个；其余保留定向检索未确认状态|
|CRC肿瘤模型基因干预报告类别|10个，包含不同强度和相反方向，不等于10个靶点|
|宿主髓系、药理救援、运输、通路内部分证据|分别1个|
|生化改变但所测迁移阴性|B4GALT2，1个|
|融合对象而非野生型|PLA2G15，1个|
|其他模型/生化或调控背景|9个|
|本轮未确认适用CRC功能|10个|
|研究工作安排|9优先核查、21保留、5暂挂|
|新增患者检验、DepMap|均未运行；不填虚假0分|

### 主要判断变化
1. **SLC6A6**：多CRC模型中有遗传干预、牛磺酸摄取、存活/耐药及肿瘤启动证据（E01）。
   这是已知机制参照，不是CAMP新靶点发现。
2. **GSTA4**：癌细胞敲除与宿主巨噬细胞/炎症研究支持不同层面的功能，需要区分细胞来源及正常炎症代价（E02–E03）。
3. **UCKL1**：原始研究明确将所报道的铁死亡机制与UMP/CMP合成作用分离（E04）。
   不能把尿苷相关性直接解释为抑瘤机制。ENO2也有非糖酵解功能背景（E09–E10）。
4. **HDC**：宿主缺失可促进炎症相关癌变，属于保护性髓系背景（E05），不能机械主张抑制。
5. **KMT5A/SETD8与AQP9**：存在模型、底物、治疗或终点不同的方向冲突（E11–E16的相应记录），应先解释而非平均成一个“支持分数”。
6. **B4GALT2**：有酶活变化，但研究中其SW480迁移试验未改变；同文其他基因的迁移结果不能转记（E23）。
7. **PLA2G15**：CRC干预对象是NFATC3–PLA2G15融合，不能作为野生型PLA2G15—GPC支持（E24）。
8. **NME6/SLC6A17**：分别有线粒体复合物/神经转运的生化证据；当前未确认可适用的CRC肿瘤功能。
   不因原相关较强就列为已证实靶点（E25–E30）。

## 新手解释
“功能有报道”表示在具体实验模型中改变了基因或相关作用，不表示该基因在所有患者都应该被抑制。
“没有确认适用文献”不表示没有功能。“优先核查”表示有具体问题值得马上解决，不是最优治疗排名。
本轮将**证据类型**与**值得继续研究的问题**分开，避免只保留P较小或文献最多的基因。
原39条q仍为不显著；文献不会改变它们的统计状态。

## 限制/反证
- 是定向证据核查，不是系统综述；未穷尽所有历史/在审文献。
- 每条记录写明“摘要”或“选读正文/图注”的实际深度，未宣称所有论文全文逐图复核。
- 部分摘要不能完整确定细胞株/独立基因终点，已明确留缺口；BCAT2的上游circRNA表型不得全部转归其单基因。
- 未进行系统撤稿数据库和图像完整性审计。GSTO2做过撤稿关键词核查，未见提示不等于认证。
- 未完成对应代谢物在患者中的因果介导验证、独立临床复现或治疗安全性评估。
- 39关系Excel字段`rho_display/p_display/q_display`来自原README四位小数摘录，不是新的高精度统计。
  高精度原值仍以固定提交TSV为准，不能用展示列覆盖原列。
- 本批未直接重新验证全部人源反应和运输底物；CANT1/ENTPD3的IMP角色、SLC29A3区室等仍待核查。
- 基因级功能注释可回接完整974行，但清楚标记哪些不是最初39条，不声称这些其他关系已被逐对功能验证。

## 当前决定
- 优先核查（9）：AQP9, BCAT2, GSTA4, HDC, KMT5A, PRMT7, SLC38A3, SLC6A6, UCKL1
- 暂挂当前代谢轴/通用治疗方向（5）：B4GALT2, BPGM, ENTPD3, PGAM2, PLA2G15
- 其余21个保留。暂挂不删除，不更改原效应和P/q。
- 05_FUNCTION与07_INTEGRATION整阶段均为PARTIAL；35基因定向核查这个批次完成。
- DepMap仍NOT_RUN/DEFERRED，所有模型数和分数留空。
- GitHub本轮创建分支实际返回403：`Resource not accessible by integration`。
  **没有远程上传成功、没有新提交**。结果在本交付包；不要把它当作main或COAD分支已有文件。

## 下一步
只补能够改变判断的那一环：模型/方向冲突、催化与非催化功能、实际分子对象、细胞来源或独立患者支持。
没有独立数据时不围绕同一批相关性反复调整到显著；候选池不缩成“只剩9个”。

## 文件与复现命令
- `COAD_35_gene_functional_review_v1.xlsx`：六张中文工作表。
- `results/gene_function_review.tsv`及JSON：35基因完整功能比较。
- `results/study_evidence.tsv`及JSON：40条来源、干预、终点、反证、阅读深度。
- `results/search_coverage.tsv`：35基因检索主题记录；不是搜索引擎完整逐请求日志。
- `source/nominal39_display_extract.tsv`：有明确精度标记的展示摘录。
- `results/nominal39_functional_display.tsv`：展示摘录与功能注释。
- `code/integrate_functional_review_nominal35_v1.py`：只读原39/974汇总表，精确核对Git blob版本，逐字段保持原统计并追加注释。
- `analysis_spec.json / source_manifest.tsv / validation.json / upload_status.json`：范围、来源与实际验证/上传状态。

在有相应COAD快照的本地仓库中执行（不需要患者矩阵）：
```sh
python code/integrate_functional_review_nominal35_v1.py --self-test
python code/integrate_functional_review_nominal35_v1.py --repo "<本地仓库>" --evidence-root "<本包解压根目录>" --out "<本地仓库>/results/COAD/05_FUNCTION/<新的唯一RUN_ID>"
```
脚本拒绝覆盖已有目录。输出原精度39与974行附加注释表，以及待审核的阶段索引追加记录；
不自动commit/push，不修改其他癌种或现用索引。
本环境已经通过脚本语法和模拟小表的字段保留测试；**真实仓库974行回接尚未执行**。
人工文献判断由版本化JSON保留，运行脚本只复现整合，不冒充自动重做文献阅读。

### 核心原文入口
- E01 SLC6A6（2014）：https://www.nature.com/articles/srep04852
- E02 GSTA4（2022）：https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2022.887127/full
- E03 GSTA4（2025）：https://www.tandfonline.com/doi/abs/10.1080/19490976.2025.2451090
- E04 UCKL1（2023）：https://pubmed.ncbi.nlm.nih.gov/37343364/
- E05 HDC（2011）：https://pubmed.ncbi.nlm.nih.gov/21170045/
- E06 NNMT（2014）：https://pubmed.ncbi.nlm.nih.gov/25201588/
- E07 NNMT（2016）：https://pubmed.ncbi.nlm.nih.gov/27323852/
- E08 NNMT（2021）：https://pubmed.ncbi.nlm.nih.gov/34539890/
- E09 ENO2（2021）：https://pubmed.ncbi.nlm.nih.gov/33934428/
- E10 ENO2（2022）：https://pubmed.ncbi.nlm.nih.gov/35954207/
- E11 KMT5A（2023）：https://pubmed.ncbi.nlm.nih.gov/36921492/
- E12 KMT5A（2026）：https://www.nature.com/articles/s41418-026-01831-5
- E13 AQP9（2017）：https://pubmed.ncbi.nlm.nih.gov/28640255/
- E14 AQP9（2023）：https://academic.oup.com/gastro/article/doi/10.1093/gastro/goad033/7205477
- E15 AQP9（2025）：https://pubmed.ncbi.nlm.nih.gov/40210882/
- E16 PRMT7（2023）：https://pubmed.ncbi.nlm.nih.gov/37541527/
- E17 SLC38A3（2024）：https://pubs.acs.org/doi/abs/10.1021/acsomega.4c00901
- E18 BCAT2（2025）：https://pubmed.ncbi.nlm.nih.gov/39642324/
- E19 ICMT（2005）：https://pubmed.ncbi.nlm.nih.gov/15784746/
- E20 ICMT（2016）：https://pubmed.ncbi.nlm.nih.gov/27710830/
- E21 GSTO2（2023）：https://onlinelibrary.wiley.com/doi/abs/10.1155/2023/4931650
- E22 SLC29A1（2015）：https://www.jstage.jst.go.jp/article/bpb/38/8/38_b14-00622/_html/-char/ja
- E23 B4GALT2（2016）：https://www.nature.com/articles/srep23642
- E24 PLA2G15（2019）：https://pubmed.ncbi.nlm.nih.gov/29909608/
- E25 NME6（2021）：https://pubmed.ncbi.nlm.nih.gov/34789336/
- E26 NME6（2023）：https://pubmed.ncbi.nlm.nih.gov/37770567/
- E27 NME6（2023）：https://pubmed.ncbi.nlm.nih.gov/37439264/
- E28 NME6（2024）：https://pubmed.ncbi.nlm.nih.gov/39273527/
- E29 SLC6A17（2008）：https://pubmed.ncbi.nlm.nih.gov/18768736/
- E30 SLC6A17（2009）：https://pubmed.ncbi.nlm.nih.gov/19147495/
- E31 HAGH（2026）：https://www.mdpi.com/2076-3921/15/2/171
- E32 NT5C3A（2014）：https://pubmed.ncbi.nlm.nih.gov/25000516/
- E33 SLC15A3（2026）：https://www.jci.org/articles/view/199709
- E34 SLC15A3（2025）：https://pubmed.ncbi.nlm.nih.gov/39719710/
- E35 SLC7A8（2018）：https://pubmed.ncbi.nlm.nih.gov/30419950/
- E36 UAP1（2015）：https://pubmed.ncbi.nlm.nih.gov/25241896/
- E37 UPP2（2011）：https://pubmed.ncbi.nlm.nih.gov/21855639/
- E38 GSTT2（2007）：https://pubmed.ncbi.nlm.nih.gov/17250773/
- E39 ENTPD3（2024）：https://pubmed.ncbi.nlm.nih.gov/38755598/
- E40 ASPA（2016）：https://pubmed.ncbi.nlm.nih.gov/27705916/
