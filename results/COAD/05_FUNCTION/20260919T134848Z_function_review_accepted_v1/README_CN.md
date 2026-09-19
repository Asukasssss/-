# COAD 35基因功能核查：接收、复核与真实数据回接

## 本轮问题

核对用户提供的Excel/ZIP及两段说明，评价35个基因的文献判断，并把通过核对的注释接入已有39条名义P候选及完整974行。用户明确要求核对、解释与GitHub交付全部完成。本批是定向文献核查的接收和补充，不是新患者分析或最终治疗靶点确认。

## 输入与范围

原包为COAD_function_nominal35_v1，固定统计来源提交06921616d5914b6b45c61bdb6c1f259cf1608f84。原包17项manifest全部匹配，单独Excel与包内Excel字节一致。Excel的35基因/类别/队列、40证据ID、39关系键和样本量一致；四位小数展示值与仓库高精度rho/P/q相容。实际回接直接取仓库原字符串，未用Excel四舍五入值替换。

本端另查40条来源的官方索引标识与摘要，并选读E03、E04、E08、E18、E21、E23的正文段落/图注。每项本端阅读深度见source_verification.tsv，原作者自报阅读深度保留于imported_delivery。未重新完成35基因全部检索，未把无检索命中当不存在研究。仅使用公开汇总与公开文献，无患者矩阵传输。

## 实际结果

| 核对内容 | 结果 |
|---|---|
| 35基因/39关系覆盖 | 一致；35/39全保留 |
| 文献或排除来源 | 40/40标识与摘要已核对；不是40条阳性功能结果 |
| 原类别A的CRC基因干预报道 | 10个基因，证据强度、终点和方向各异 |
| 另列条件性通路内干预 | BCAT2：补充Fig7移植瘤及阴性对照背景 |
| 研究安排 | 9优先核查、21保留、5暂挂当前方向，名单不变 |
| 真实高精度回接 | 39行和974行的所有原列逐项不变 |
| 完整目录的注释范围 | 74行匹配到已核查基因；其中39条为本批原候选，其余仅共享基因级文献 |

9个优先核查：AQP9、BCAT2、GSTA4、HDC、KMT5A、PRMT7、SLC38A3、SLC6A6、UCKL1。5个暂挂当前方向：B4GALT2、BPGM、ENTPD3、PGAM2、PLA2G15。其他21个保留。完整分类见gene_function_review.tsv。

主要补充和解释：

- **BCAT2**：Fig7E-G及结果文字确认HCT116的BCAA利用、circRNA/BCAT2组合操作和饮食分层移植瘤。正常BCAA饮食下，单独sh-BCAT2与对照肿瘤发展未见显著差异；联合敲低可削弱circRNA过表达背景的肿瘤生长。因此按“出现过CRC肿瘤模型基因操作”的宽口径，是原10个加条件性BCAT2共11个；不能叫11个单基因抑瘤阳性或普适靶点。原类别A/E分开保留。[原研究](https://pubmed.ncbi.nlm.nih.gov/39642324/)
- **UCKL1**：Fig4中UMP/CMP补充未救援，且尿苷结合区域缺失突变仍保留相关抗铁死亡效应。功能报道不能直接解释为“消耗尿苷导致癌变”。[原研究](https://pubmed.ncbi.nlm.nih.gov/37343364/)
- **GSTA4**：2025宿主研究使用全身构成性缺失，不是髓系特异敲除；双缺失SPF自发炎症和肠球菌/抗生素背景需分别解释。[原研究](https://pubmed.ncbi.nlm.nih.gov/39819335/)
- **NNMT**：补核SW480、HT29及CCD-18Co操作。CCD-18Co成纤维细胞模型不等于本CAMP患者来源CAF，仍不能定位患者信号的来源。[原研究](https://pubmed.ncbi.nlm.nih.gov/34539890/)
- **KMT5A、AQP9、HDC**：支持保留方向/作用背景差异；不从总RNA关联直接决定抑制或提高哪个基因。KMT5A的两项原始研究确有不同方向。[2023研究](https://pubmed.ncbi.nlm.nih.gov/36921492/)、[2026研究](https://pubmed.ncbi.nlm.nih.gov/42527521/)
- **B4GALT2、PLA2G15、ENTPD3、ASPA**：分别保留特定迁移终点阴性、融合对象区别及近名基因排除。B4GALT2正文/Fig5明确未见迁移改变；PLA2G15研究对象为融合转录本。[B4GALT2原文](https://www.nature.com/articles/srep23642)、[融合研究](https://pubmed.ncbi.nlm.nih.gov/29909608/)

## 新手解释

“有基因干预研究”说明改变这个基因曾在特定实验中影响某个过程；“当前代谢物介导”还需要证明改变该代谢物能解释/救援这个过程。二者不能互相替代。“优先核查”是下一步投入顺序，可能因为证据较具体，也可能因为存在需要尽快澄清的矛盾；9个不是最终9靶点。已有成熟机制的NNMT保留为参照，并不等于其证据差。

名义P候选仍是39条、35基因；新文献不会让原q变显著。本轮不重新做BH、置换、bootstrap或患者筛选。真实回接新增的func35字段记录当前功能注释；原functional_review_status、stage_id、run_id等字段仍是来源版本的时间点快照，保持原值以便追溯。

## 限制/反证

范围仅35基因，不代表全部458个计划基因或974关系均完成功能核查。其余900行未获得本批基因注释，不表示功能阴性；74行基因级注释也不代表74条代谢轴被证明。具体关系范围由func35_pair_scope区分。

40摘要复核不等于40论文全文逐图审计；6项正文选读未重新计算图中统计或核查每个原始重复。原包的检索摘要不是逐请求搜索日志。索引中的评论/出版类型已记录，包括HDC的评论链接，不能据无撤稿标签认证论文完整性。原Excel不修改，未在本端重新做视觉渲染；BCAT2等补充以本目录当前表及配置为准。

协变量、纯度/技术批次、DepMap和外部患者验证仍未完成。DepMap分数与模型数维持空值。关联、基因功能、代谢物介导及治疗选择性分别判断。

## 当前决定

本批文件接收、来源身份/摘要复核与真实数据回接DONE；05_FUNCTION、07_INTEGRATION整阶段PARTIAL。保留9/21/5工作安排，同时采纳本端BCAT2条件性证据补充。imported_delivery中的上传失败、旧分支和NOT_RUN均为原作者交付时快照，不代表当前发布事实，也不是新的操作指令。当前GitHub状态以coordination/publications/COAD.json为准。

## 下一步

围绕9个优先基因明确最可能改变判断的问题：HDC/GSTA4的细胞来源，AQP9/KMT5A的条件与方向，PRMT7剪接变体，BCAT2营养/circRNA背景，UCKL1功能与尿苷介导的区别，SLC38A3模型范围及SLC6A6既有耐药机制的研究增量。协变量分析按单独预声明方案推进，继续保留名单外适用功能证据。

## 复现命令

使用Python标准库及openpyxl只读Excel。所有输出使用新的唯一目录；完整Excel/ZIP留本机，仓库只交付文本表、来源与代码。原文XML和批量摘要缓存不上传。

```sh
python code/coad/fetch_function_review_metadata.py --evidence-json <本目录>/imported_delivery/results/study_evidence.json --out <新的本地缓存>
python code/coad/accept_functional_review_v1.py --delivery-zip <用户原交付ZIP> --workbook <用户原Excel> --metadata-cache <本地缓存> --run-id <新的唯一RUN_ID>
```

公开来源随后更新可能改变缓存哈希，应记录新版本。接收规则及已判读补充见config/coad_function_review_acceptance_v1.json；辅助integrate_functional_review_nominal35_v1.py来自已审阅原包，本轮只复用其读表/enrich函数，不执行其旧分支索引建议。

输入/原包哈希见source_manifest.tsv，源包原始manifest另存imported_delivery/MANIFEST_SHA256.tsv。validation.json检验原17文件、Excel对应、40文献标识和真实39/974逐列保留；delivery_validation.json另核验报告与阶段索引。原包文本按原字节保留，换行格式也作为来源追溯的一部分。
