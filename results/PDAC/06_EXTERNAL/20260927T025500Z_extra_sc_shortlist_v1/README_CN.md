# PDAC 新增单细胞数据初筛（2026-09-27）

## 本轮问题
为LYPLA1、LYPLA2及完整候选池寻找现有三套之外的补充来源数据。筛选基于设计、数据开放性和注释潜力，不看目标基因是否显著。

## 输入与范围
核对GEO官方索引、原始研究及作者资源。现有已分析队列为GSE263733、GSE278688、GSE242230。以下是研究设计层面的优先级，不是完成矩阵QC后的质量认证。

## 实际结果

| 队列 | 公开设计 | 适合补充的问题 | 限制 |
|---|---|---|---|
| GSE155698 | 16份PDAC组织、3份邻近正常；另有PBMC | 全细胞来源及邻近正常上皮对照 | 11A/11B要解析；不能把17个肿瘤文库当17患者；正常少、不假定配对 |
| GSE202051/SCP1089 | 43份原发PDAC：18未治疗、25新辅助治疗；另2非恶性标本 | 恶性上皮来源；优先未治疗部分 | 单核；74条GEO记录含重复等，不等于74供者；两个入口不是两队列 |
| GSE217845 | 10个肿瘤文库、6个外周血文库 | 髓系/TAM内部来源 | 至少1个肿瘤文库明确髓系富集；不能直接合并算组成；无正常胰腺设计 |
| GSE291124 | 作者报告17名未治疗原发PDAC患者 | 较新的未治疗补充队列 | 单核；无正常胰腺；细胞标签仍需连接 |
| GSE154778 | 10份原发、6份转移活检 | 备用来源队列 | 无正常组织；10x版本不同；原发和转移分开 |

推荐先接GSE155698和GSE202051的未治疗部分，针对髓系再接GSE217845；GSE291124作追加，GSE154778作备用。此排序是本轮针对研究问题的判断，不是按期刊或细胞数排名。

## 新手解释
公开矩阵、论文中的细胞类型图和能与barcode连接的作者标签是三件事。GSE155698有作者代码，GSE202051有H5AD及SCP门户；但本轮尚未验证逐细胞标签及计数层。GSE217845的TAM研究定位适合补充亚群，并不意味着有“恶性髓系”标签。邻近正常组织也不能自动当健康胰腺或与任意肿瘤配对。

## 限制/反证
部分GEO/Nature直连遇到CAPTCHA/跳转，使用官方检索索引与原始研究交叉核对。SCP1089页面实际显示88,031个细胞；该数字仅对应这个门户，不推广到整个GSE202051研究。服务器SSH连接超时，未检视服务器矩阵，未下载原始数据，未计算任何新增LYPLA1/LYPLA2结果。完整患者重叠检查、每患者细胞数、双细胞、基因覆盖和注释连接均未完成。GSE155698/GSE154778/GSE202051在旧选择表中已是备选或补充，不冒称从未发现；本轮新增重点为GSE217845与GSE291124。

同一研究的GEO、SCP及整合图谱不重复计为验证；GSE217845是GSE217847的人类scRNA子系列。GSE291124原始17人数据与其论文重用的Steele队列分开。不同研究不等于患者独立性已逐例证实。

## 当前决定
本轮公开资料初筛DONE；矩阵QC、逐细胞标签连接及新增基因分析NOT_RUN；独立性NEEDS_REVIEW。三优先、一补充、一备用。不能保证增加队列后会获得显著结果。

## 下一步
实际接入时先核验作者单位、原发/转移、治疗和富集策略，再按作者单位汇总完整候选池。正常上皮及髓系组须报告实际合格供者数，不用总细胞量替代独立单位。

## 复现命令
python code/pdac/extra_sc_shortlist_v1.py

脚本只重建人工核对的筛选记录，不自动抓取数据或重现论文分析。完整入口、设计、下载格式、注释缺项见cohort_shortlist.tsv。source_manifest.tsv中的远程网页SHA256为NA，因为未保留网页字节；不编造哈希。

## 原始来源

- [GSE155698](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155698)；[原研究或作者资源](https://github.com/PascaDiMagliano-Lab/MultimodalMappingPDA-scRNASeq)。

- [GSE202051/SCP1089](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202051)；[原研究或作者资源](https://singlecell.broadinstitute.org/single_cell/study/SCP1089/human-treatment-naive-pdac-snuc-seq)。

- [GSE217845](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217845)；[原研究或作者资源](https://pubmed.ncbi.nlm.nih.gov/37914939/)。

- [GSE291124](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291124)；[原研究或作者资源](https://pmc.ncbi.nlm.nih.gov/articles/PMC12866174/)。

- [GSE154778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154778)；[原研究或作者资源](https://pmc.ncbi.nlm.nih.gov/articles/PMC7523332/)。
