# BRCA：LYPLA1高低表达癌上皮的差异基因与富集

## 本轮问题
在恶性上皮中，LYPLA1高表达细胞伴随哪些基因和转录通路变化？按用户指定以细胞直接比较，不做供者汇总；保留供者/亚型/深度构成检查。不是LYPLA1干预实验。

## 输入与范围
Wu作者恶性上皮注释与Pal数据的Chen2026重注释分别分析；Pal仅batch=pal_2021。每套数据先将LYPLA1按细胞总UMI标准化为log1p(CP10k)，在已检出细胞内固定中位数。主比较：高于中位数 vs 大于0且不高于中位数；补充：高表达 vs 未检出。相同值不拆分，无阈值搜索，无供者细胞数门槛。

### 分组
|cohort|group|n_cells|positive_median_cut|n_donors|median_counts|median_detected_genes|
|---|---|---|---|---|---|---|
|Wu2021|high|7222|1.064761516|20|9274.5|2586.5|
|Wu2021|low_detected|7222|1.064761516|20|17516.5|3729.0|
|Wu2021|zero|10045|1.064761516|20|4583.0|1633.0|
|Pal2021_reprocessed|high|25622|1.21676196|33|4449.5|1641.0|
|Pal2021_reprocessed|low_detected|25623|1.21676196|31|10333.0|2716.0|
|Pal2021_reprocessed|zero|78833|1.21676196|33|2274.0|947.0|

## 实际结果
### 全量差异与富集覆盖
|cohort|contrast|n_high|n_reference|tested_genes|P_lt05|q_contrast_lt05|up_P_effect|down_P_effect|GSEA_terms|GSEA_nominalP_lt05|GSEA_empiricalFDR_lt025|
|---|---|---|---|---|---|---|---|---|---|---|---|
|Wu2021|high_vs_low_detected|7222|7222|14558|10910|10713|1358|1481|1037|329|445|
|Wu2021|high_vs_zero|7222|10045|14558|10361|10070|1305|1568|1037|332|378|
|Pal2021_reprocessed|high_vs_low_detected|25622|25623|13204|11235|11170|1191|869|1004|358|482|
|Pal2021_reprocessed|high_vs_zero|25622|78833|13204|12462|12454|1491|282|1004|416|551|

已保存4份全基因差异表、4份Reactome GSEA、4份GO/Reactome上下调ORA及跨队列对照。GSEA每项保留leading-edge基因；ORA保留命中基因，绝不只交付漂亮通路。

主比较两研究均见翻译/核糖体相关表达程序负向富集。它们是多个高度重叠的条目，不能算多次独立发现，也不能直接解释为翻译速率降低。

Wu中“Peroxisomal Lipid Metabolism”正向富集（NES约1.811，经验FDR约0.055）；1000次置换的名义P在软件中报告为0，意为本次置换未出现更极端结果，绝非真实P=0。Pal同条目NES约0.931，P约0.607，未支持。因此不能将Wu脂质线索称为两队列一致。

Wu该条目的leading-edge包括CROT、SLC27A2、CRAT、HSD17B4、ACOT4、ACOX2、ACBD4、ACOX1、ALDH3A2、SYCP2、SLC25A17、ACOX3；成员命中不意味着全部是直接代谢酶或LYPLA1下游靶点。完整精确数值见GSEA表。

### 深度背景
|cohort|contrast|metric|mean_high|mean_reference|standardized_mean_difference|
|---|---|---|---|---|---|
|Wu2021|high_vs_low_detected|log1p_total_counts|9.1261|9.7783|-0.7785|
|Wu2021|high_vs_low_detected|detected_genes|2894.505|3886.4585|-0.6722|
|Wu2021|high_vs_zero|log1p_total_counts|9.1261|8.4093|0.7876|
|Wu2021|high_vs_zero|detected_genes|2894.505|1850.3353|0.785|
|Pal2021_reprocessed|high_vs_low_detected|log1p_total_counts|8.4442|9.2934|-1.1302|
|Pal2021_reprocessed|high_vs_low_detected|detected_genes|1948.6433|2921.7786|-0.8558|
|Pal2021_reprocessed|high_vs_zero|log1p_total_counts|8.4442|7.7112|0.7725|
|Pal2021_reprocessed|high_vs_zero|detected_genes|1948.6433|1163.441|0.7739|

两研究高组总UMI及检出基因数的中位数均低于低组。因此高LYPLA1归一化值可能夹杂库大小、复杂度、细胞状态和患者构成影响。主比较排除0值不能自动消除这些影响。Wu高组ER+细胞占70.49%，低组占39.05%；高组TNBC占22.51%，低组占51.12%。因此ESR1和雌激素表达信号可能受到亚型混杂，不能当作LYPLA1独立效应。亚型、治疗构成见subtype_treatment_composition.tsv；Pal缺相应配置，未猜测补充。供者身份留服务器，仅公开人数、占比摘要和构成差异指标。

## 新手解释
这些结果描述“LYPLA1高的癌细胞还表现出什么特征”。当前最一致的是部分翻译相关基因偏低；脂质相关正向线索主要见于Wu，不应一概写成脂质代谢增强，更不能写成LYPLA1驱动某条通路。

## 限制与反证
1. 同患者细胞不独立，本轮按用户要求不调整供者，细胞P可能过小，不能替代患者层面显著性。高/低由表达定义，也不是随机分组。
2. 测序深度、患者和亚型构成仅检查，未调整；当前不能将差异全部归因于LYPLA1。
3. GSEA使用基因集置换，不是供者/细胞标签置换；NES和经验FDR不是独立患者验证。零P/零FDR只保留软件原输出并说明置换分辨率。
4. 未检出不代表不存在；高vs零与高vs低是不同问题。两研究数据处理及注释不同，交集结果仍为描述性一致性。
5. Reactome_2022、GO_Biological_Process_2023为固定版本，非声称最新。固定15至500个已检验基因的通路才能评估；小集合被过滤不等于阴性，见GSEA_pathway_coverage.tsv。
6. 未开展新模型调参、再分亚型检验、通量推断或机制实验。旧结果不覆盖。

## 当前决定
保留全部差异基因和通路作为LYPLA1相关状态的探索结果；Wu脂质条目可作为限定线索，不作为跨队列稳定机制。优先理解深度与亚型构成，再决定是否值得独立验证；不继续换阈值争取显著。

## 下一步
本轮交付到“基因与通路名单、驱动基因、两研究对照和限制”停止。已有CAMP LYPLA1—LPC直接关系单独保留，不因本轮富集支持有限而删除，也不把它与富集拼成已证实因果链。

## 复现命令与字段
先在server165建新独占运行目录及.running，将三个分析/展示/审计脚本复制进去。
`python3 brca_lypla1_highlow_v1.py <RUN>`
`python3 brca_lypla1_highlow_audit_v1.py <RUN>`
`python3 brca_lypla1_highlow_present_v1.py <RUN>`
`python3 brca_lypla1_highlow_report_v1.py <RUN>/public`

方法：稀疏计数按细胞CP10k后log1p；全癌上皮检出率>=1%、唯一符号、对比有非零方差。LYPLA1自身不测试、不富集。双侧Welch t检验用于细胞探索，提供均值差、检出率和描述性log2((mean CP10k高+0.1)/(mean CP10k参照+0.1))，后者带平滑，不能当作原始精确倍数。GSEA使用全部可检验基因按带符号t排序，1000次置换，seed20260925，weight1，Reactome。ORA用P<0.05且描述性log2比值绝对值>=0.25的上下调集合，背景为该比较全部已检验基因，超几何右尾，包含零命中条目。阈值是候选展示规则，不是效应超过阈值的正式检验。

DE的q_contrast=本比较全基因BH；q_joint4=四比较全基因联合BH。ORA的q_contrast=本比较两个方向×两个库全部条目BH；q_joint4=四比较联合BH。GSEA经验FDR由包计算，不能与BH列混为一谈。重点展示P但全q留档。

明细表字段映射：cohort、contrast、gene为统计键；n_high为n，n_reference为参照细胞数；mean_log_difference和descriptive_log2ratio是两种独立效应尺度；p_value是细胞P。源路径、SHA和数据库URL见source_manifest.tsv。逐细胞信息只在服务器private中保存。审计独立复算部分基因Welch统计、部分GSEA富集分数及ORA精确概率；没有声称独立重跑全部置换。图表只展示固定规则选取的部分，完整表是解释依据。

GSEA实现参考：https://gseapy.readthedocs.io/en/latest/run.html 。本运行实际版本及参数以analysis_spec.json为准。
