# BRCA：Tang外部患者关联、Oslo2可评估性与117基因RNA背景

## 本轮问题
让FUSCC之外的配套资料评估原174条关系，并将RNA单层背景与代谢物—RNA关联分开。完整117候选保留，不以新队列是否显著设入场门槛。

## 输入与范围
Tang2014 Table S2实际包含399个代谢特征、25份肿瘤和5份正常。逐病例代谢值沿用作者归一化和填补结果，没有重新填补。以作者TCGA Designation字段去除首尾空格、补标准TCGA-前缀，与cBioPortal brca_tcga_pub2015的patientId精确连接，只纳入Primary Solid Tumor且每病例唯一标本。所选RNA profile为brca_tcga_pub2015_rna_seq_v2_mrna（RSEM连续值，非z-score）。当前可连接20例，不等同于论文当年23例；不是全TCGA具有代谢组，也不宣称对应同一分装。该RNA版本不冒充原文2014统计版本。

Oslo2 Additional file2的Add. Table2已实际读取：228行，18项NMR代谢物。GSE58212的283份样本元数据使用OSL2编号，代谢物表使用MCOS编号。GEO完整元数据与后续2017研究附件12张表均未发现MCOS对照。没有按排序、数字后缀、分型或表达模式猜配对，因此当前连接数0，RNA完整矩阵暂不追加下载。

GSE42568通过官方GEO逐样本处理表取得121份记录（104癌、17正常），复用作者log2 GC-RMA。GPL570唯一Entrez注释对应的全部探针对每个样本取中位数；不按显著性或相关强弱挑探针，不重新标准化。每组作为独立作者标本处理，未推断配对，不能认为新增核实了所有患者身份及批次。

## 实际结果
- Tang：73个代谢物键中66个唯一作者名称对应；1项乳酸立体范围待核查，6项未有唯一名称匹配。原174关系全部保留，163条实际可算，均n=20；2条BH q<0.05。
- Oslo2：8个代谢物键可建立明确名称对应、另2个需核查化学范围；对应34条关系仅受阻于样本对照，另4条还有化学范围限制；136条不在18项面板内。所有174条仍有去向。
- GSE42568：117候选全部保留，112个可评估，57个RNA差异q<0.05。AMDHD2、CPT1B、PAH、PRODH、PYCR2无满足当前规则的唯一Entrez探针，保留不可评估，不作表达阴性。

|Tang关系|ρ|点95%区间|本轮q|
|---|---:|---|---:|
|MDH1—苹果酸|+0.746|+0.441～+0.920|0.0163|
|PNP—次黄嘌呤|+0.687|+0.361～+0.859|0.0489|
|GPI—G6P|−0.129|−0.559～+0.372|0.9151|
|GPI—F6P|−0.203|−0.663～+0.288|0.7790|
|GPCPD1—GPC|−0.328|−0.676～+0.127|0.5466|
|ASNS—谷氨酰胺|−0.149|−0.611～+0.319|0.8842|
|GLS—谷氨酰胺|−0.344|−0.719～+0.108|0.5383|
|SLC6A8—肌酸|+0.382|−0.049～+0.685|0.5208|

MDH1—苹果酸在CAMP原未调相关为+0.236、FUSCC为+0.016；PNP—次黄嘌呤分别为+0.272、+0.118。它们在Tang有新支持，但不能称三个队列均显著验证；PNP—次黄嘌呤不等于之前PNP—鸟嘌呤。所有174条跨队列明细见all174_cross_cohort_context.tsv，不只展示显著关系。

## 新手解释
ρ为患者之间的共同变化程度，不是浓度倍数。GPI两条关系在这20例中的点估计为负，区间跨零，没有支持CAMP的正关联；这也是有用的结果，不能删掉或称外部复现成功。ASNS和GLS仍为负向，区间宽且本轮未显著，也不能据此推翻此前FUSCC结果。

RNA差异是“癌组织与正常组织的表达是否不同”，并非“该基因与代谢物是否有关”。例如GSE42568中ASNS的RNA差异为+0.835 log2单位、GLS为+0.046；这与患者内部谷氨酰胺关联是两个不同问题，不能互相替代。

## 限制/反证
Tang每项只有20例，区间为逐项bootstrap区间而非同时覆盖区间；原代谢组已经填补，无法恢复原缺失模式。RNA版本、分型、组织组成及分装差异仍存在；本轮未追加ER、纯度等模型，也未做癌种/队列间效应差异检验。

名称对应是作者化合物标签层面的匹配，没有重新核验谱图、标准品、异构体或反应通量。G6P使用明确别名glucose-6-phosphate (G6P)，没有以葡萄糖替代；F6P为独立作者特征。乳酸和Oslo谷胱甘肽范围歧义未强行通过。

Tang的P为9999次双侧置换加1，最小可报告P=0.0001；PNP q=0.0489接近阈值，应考虑蒙特卡洛误差，不夸大其与0.05两侧结果的区别。没有为改善显著性改检验次数或校正范围。bootstrap2000次，seed=20260921+原关系行序号。每队列原174条关系中全部可评估项一次BH（Tang163项），缺项保留但不伪造P=1。GSE42568另为112项BH；效应为未配对Welch均值差及Welch95%区间。完整117RNA预定范围不改为显著子集。

独立脚本逐项复核Tang rho、GSE42568效应/P及两套BH，数值误差均<1e-10；未重新生成置换P与bootstrap区间。历史CAMP/FUSCC/Asns统计未覆盖。源矩阵与完整样本记录只在server165。

## 当前决定
保持两条Tang新支持关系为候选证据，保留GPI未获支持及其他关系的宽区间。Oslo2目前状态为NEEDS_REVIEW，34条已匹配关系只等作者MCOS—OSL2编号对照，不能将论文201例当作已连接数。原117比较表的215个历史列完整保留，追加9列本轮结果。

## 下一步
先解决Oslo2真实编号对照（需要作者公开映射或可追溯的补充文件）；未取得前不算关联。新支持的MDH1与PNP需结合完整旧功能证据判断，不以本轮两个q替换全117候选池。Tyanova蛋白与More代谢物单层资料本轮未下载分析，作为后续可选补充，不登记已验证。

## 复现命令
运行目录：/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A/20260921T072503Z_oslo_tang_external_v1

在该服务器目录执行（source/relations174.tsv来自此前全174冻结关系；输入缓存哈希见source_manifest.tsv）：
```sh
python3 brca_external_fetch_more_v1.py "$PWD" tcga
python3 brca_external_fetch_more_v1.py "$PWD" geo
python3 brca_oslo_tang_analysis_v1.py "$PWD"
python3 brca_gse42568_background_v1.py "$PWD"
python3 brca_validate_new_external_v1.py "$PWD"
```
论文xlsx可用brca_oslo_tang_download_v1.py下载。首次GEO FTP矩阵请求超过600秒未完成，已停止自己的下载进程并改用官方GEO处理表；下载脚本已增加子进程墙钟超时。source缓存只读复用；没有下载原始CEL或质谱。

取回public汇总后，在仓库运行：
```sh
python code/brca_finalize_new_external_v1.py
python tools/check_repository.py
```

输出：Tang_all174_associations.tsv、Oslo_all174_evaluability.tsv、GSE42568_all117_RNA.tsv为统一字段顺序；all117_comparison_external_appended.tsv为保留历史的比较表；Tang_selected_relationships.png为八条关系的区间图，展示范围不改变全163项校正。

来源：
- Tang2014：https://link.springer.com/article/10.1186/s13058-014-0415-9
- Oslo2：https://link.springer.com/article/10.1186/s40170-016-0152-x
- Oslo2后续附件核查：https://link.springer.com/article/10.1186/s13058-017-0812-y
- GEO：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42568 ，GSE58212同入口。
- TCGA RNA：https://www.cbioportal.org/api/studies/brca_tcga_pub2015/molecular-profiles
