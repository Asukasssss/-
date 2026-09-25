# GBM 当前结果：CARE 2025 单细胞来源替换

## 本轮问题
按用户要求用质量更好、细胞量足够且较新的数据替换 Darmanis2017 和 Neftel2019。新主来源为 CARE 2025，旧版完整保留。

## 输入与范围
CARE 2025（Nature Genetics；GSE274546），10x 单核 RNA 测序，7家医院、59位IDHwt患者、121标本；作者质控后429,305核。细胞类型来自作者2025-01-08注释，原计数矩阵逐条码连接。采用作者患者ID，不根据T1/T2编号推断初发或配对。

- 主分析：明确未接受术前放疗/烷化药的初发样本，55患者、55标本、205,880核。
- 排除：1个术前已接受放疗与TMZ的作者初发标本，795核；规则在观察候选表达前确定。激素使用保留，未将此队列称作未接受任何药物。
- 复发补充：59患者、65标本、222,630核，多标本在同一患者内合并。
- 初发和复发不是两个独立队列；本次没有新增独立复现结论。

## 实际结果
全部142个候选均保留；精确名称覆盖136/142。未匹配：IARS1, G6PC1, GBA1, GBA3, SLC60A2, LARS1，标记不可评估，不补零、不猜别名。

LYPLA1 初发最高表达类别：**Myeloid**，次高：Vascular；整患者bootstrap保持率99.8%（有效1000/1000）。该最高类别中可评估患者51位，患者等权检出比例38.6%。复发最高类别Myeloid，保持率100.0%；属于同队列补充。

此前代谢物/bulk结果原样保留并逐列精确核查：74肿瘤、6非配对正常；697可比较代谢特征，480个P<0.05、460个当前族q<0.05；当前映射23特征→171关系→142基因。168关系可评估，14条P<0.05，0条q<0.05。RNA139基因可评估，75个q<0.05。LYPLA1的RNA差异和代谢关联P/q没有因单细胞替换而改变。

## 新手解释
细胞量增加改善细胞类别覆盖，55位主分析患者用于跨患者稳定性。表达按完整33,538基因库归一化到CP10k后log1p，先在患者×类别内平均，再患者等权；每患者类别至少20核、每类别至少3患者。来源排名至少2可评估类别，最高平均检出率至少1%。1000次整患者bootstrap，重复抽中同一患者只增加权重。

灰色/NA表示不可评估。作者Other类别保留表达描述但不参与明确来源排名。单核和完整单细胞技术、组织采样及细胞类型覆盖不同，旧/新最高类别改变不能解释为生物学变化。

## 限制/反证
CARE为可再次手术的纵向病例队列，有临床选择偏倚。细胞来源只回答转录本主要在哪些类别表达，不证明代谢物来源、酶活、代谢通量或机制。没有新增聚类、CellChat、拟时序或功能筛选。原bulk比较仍有GTEx来源混杂；457工作特征待补精确生化映射。

## 当前决定
新主来源改为CARE未接受术前放疗/烷化药的初发队列，复发单列补充。LYPLA1及全部142候选保留；当前仍没有代谢物—RNA关系通过FDR，不能因单细胞背景升级成机制证据。

## 下一步
按新的初发表达来源讨论候选及适合的细胞模型。若要声称跨队列独立重复，须另增加经患者来源去重的队列；本次未作此声明。

## 复现命令
代码：gbm_download_care_v2.py → gbm_care_aggregate_v2.R（4个互斥文件分片）→ gbm_care_profile_v2.py → gbm_validate_care_v2.py → gbm_integrate_care_v2.py → gbm_document_care_v2.py。
服务器运行：/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/GBM/A/20260925T151224Z_sc_replacement_v2。
公开路径：results/GBM/06_EXTERNAL/20260925T151224Z_sc_replacement_v2 与 results/GBM/07_INTEGRATION/20260925T151224Z_sc_replacement_v2。原计数、逐细胞与患者级资料仅在server165；当前交付为聚合结果。

## 文件与来源
- GBM_CARE2025_全候选142基因171关系.xlsx：全候选工作簿。
- figures/GBM_all142_CARE_cell_sources.pdf：全142点图；初发/复发分别标注。
- figures/GBM_CARE_primary_all142_heatmap.pdf：初发全候选相对表达热图。
- figures/LYPLA1_CARE_sources.png：LYPLA1两分区表达来源。
- [原始研究](https://www.nature.com/articles/s41588-025-02168-4)；[GEO计数](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274546)；[作者固定版本注释](https://github.com/dravishays/GBM-CARE-WT/tree/383b0a320f0993a4529167cdb591dac63e103674)。
