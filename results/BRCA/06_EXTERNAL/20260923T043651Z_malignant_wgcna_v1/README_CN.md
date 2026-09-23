# BRCA恶性上皮共表达：LYPLA1模块 v1

交付状态：DONE：100次置换保留已完成。

## 本轮问题

在作者标注的恶性上皮内构建共表达模块，查明LYPLA1的模块归属、同模块基因、功能富集和跨研究表现。研究目的为潜在靶点的数据挖掘，不是因果或治疗依赖验证。

## 输入与范围

- 发现：Wu2021原作者celltype_major；外部：Pal2021_reprocessed，使用Chen2026对Pal计数的重注释，不冒充Pal原作者标签。均复用冻结标签，未重聚类或重新判定恶性。
- 实际覆盖见coverage.tsv；Wu纳入20个作者供者标签、160个metacells；Pal纳入32个作者供者标签、256个metacells。
- 每位纳入供者至少120个恶性上皮细胞；先最多抽1000细胞用于候选矩阵/PCA，再在供者内部形成8个不重叠的15细胞近邻组。每位供者对网络贡献相同，未筛选LYPLA1高表达细胞。
- Wu按唯一基因符号、非线粒体、候选细胞检出率≥5%、至少3位供者检出、正方差过滤，选log1p归一化方差最高6000基因。LYPLA1自行通过过滤，没有强制入网；不是只用原候选基因建网。
- 候选细胞集每供者上限相同，但实际数可不同，因而基因筛选/PCA不是严格供者等权。最终构网8组/供者等权。供者汇总使用候选细胞集，最多1000/供者。
- 源矩阵、细胞成员、供者表达均留server165。source_manifest.tsv记录源SHA256。原CAMP统计不变。

## 实际结果

|项目|结果|
|---|---|
|构网基因|6000|
|有效模块（不含grey）|17|
|软阈值|12，fallback=False|
|LYPLA1模块|pink|
|模块基因数|213|
|LYPLA1 Wu kME|0.4974|
|模块内kME名次|127.0|
|Pal固定模块LYPLA1 kME|0.4982|
|留出一位供者的kME范围（模块排除LYPLA1自身）|0.401～0.555|
|供者层面外部Zsummary|1.321|
|本模块GO BP联合q＜0.05条目|48|

完整模块和全部基因均交付，不只展示LYPLA1所在模块。模块名称仅为颜色标签。

GO前10项（未显著仍展示）：

|GO生物过程|重叠基因数|全部模块×条目联合q|
|---|---|---|
|purine ribonucleoside triphosphate biosynthetic process|15|7.39e-05|
|proton motive force-driven ATP synthesis|13|7.46e-05|
|purine nucleoside triphosphate biosynthetic process|15|8.58e-05|
|ribonucleoside triphosphate biosynthetic process|15|9.9e-05|
|ATP biosynthetic process|14|0.000121|
|proton motive force-driven mitochondrial ATP synthesis|12|0.000132|
|nucleoside triphosphate biosynthetic process|15|0.000235|
|ribonucleotide biosynthetic process|18|0.000253|
|purine-containing compound biosynthetic process|19|0.000303|
|ribose phosphate biosynthetic process|18|0.000445|

脂质关键词索引是方便查阅的附加索引，不改变检验范围或q；只包含名称匹配，不是穷尽的脂质功能判定：

- cellular lipid catabolic process：重叠8基因，联合q=0.337。
- fatty acid beta-oxidation：重叠5基因，联合q=0.378。
- lipid catabolic process：重叠9基因，联合q=0.42。
- fatty acid oxidation：重叠5基因，联合q=0.668。
- lipid oxidation：重叠5基因，联合q=0.698。
- fatty acid catabolic process：重叠5基因，联合q=0.762。
- membrane lipid catabolic process：重叠2基因，联合q=1。
- lipid modification：重叠5基因，联合q=1。

## 新手解释

模块表示RNA经常一起变化的一组基因。kME表示某个基因与该模块概括表达的相关程度；不是表达倍数，也不是干预效应。中心性高只能提示值得进一步比较，不能因此称为驱动靶点。

GO富集表示该模块相较本次构网背景含有更多某类注释基因；不等于在细胞里测到相应代谢物、酶活或通量。若脂质条目没有通过校正，不能因为LYPLA1已知生化功能就把整个模块命名为脂质模块。

外部Zsummary用于探索性共表达结构保留，常用2和10参考线，但它受模块大小、样本覆盖影响，不能视为P值或治疗验证。Pal kME只是固定模块成员的描述性相关，不是重新发现了相同模块。

## 限制与反证

- 使用自定义不重叠metacells加WGCNA 1.74，未运行hdWGCNA软件包。邻居在同供者内选，PCA分别按研究计算；不存在跨研究混合metacell。
- bicor遇到零MAD的列按预设pearsonFallback=individual回退至Pearson，软件日志中的相应警告被保留。
- 160/256个细胞组不是160/256名独立患者。主网络用于发现；正式外部保留以供者汇总作为输入（Wu 20、Pal 32）。作者身份标签未做基因型独立性确认。
- 留出一位供者只重算固定模块的概括表达与LYPLA1相关，没有重建全网络，所以不能称为模块归属重采样稳定率。
- 没有消除亚型、供者、技术批次、细胞周期或CNV的全部影响。不能把跨细胞共表达直接解释为细胞内调控；亚型均值表为描述性，无新亚型P/q。
- Pal无唯一同名记录的基因单列，不拼接别名、不填0；供者保留进一步排除零方差基因。模块保留100次置换，大模块最多抽1000基因，属于探索性精度。
- 无正常上皮参与建网，本轮不重新检验恶性比正常是否高表达。
- GO版本与全部检验族见enrichment_spec.json；BH跨全部模块×全部合格BP条目（包括零重叠）。关键词索引不是单独校正的脂质检验族。
- 已知LYPLA1乳腺癌研究及原患者结果仍是独立证据栏，本轮不改写已有研究的新颖性，也不把共表达补充当作功能实验。

## 当前决定

保留LYPLA1作为已有患者/表达线索上的共表达专题对象。是否与脂质程序相关，以完整GO结果、同模块基因和外部保留共同解释；不为得到某个模块反复改阈值。其他模块及未支持条目完整保留。后续无需自动展开完整机制研究。

## 下一步

先阅读本模块前列基因、GO全表与四张图，判断这条模块信息是否真的增加候选价值。若主要反映普遍翻译/增殖程序，也如实记录，不强行包装成脂质特异性。后续干预或其他队列是独立新问题，不是本轮完成条件。

## 复现命令与记录

在server165新建独占运行目录，复制脚本、冻结的sc117读取helper和Wu/Pal配置；安装WGCNA到运行目录Rlib。配置文件及helper均在仓库code/。依次运行：

```text
Rscript brca_malignant_wgcna_install_v1.R RUN
python3 brca_malignant_metacells_v1.py RUN
Rscript brca_malignant_wgcna_v1.R RUN
Rscript brca_malignant_wgcna_followup_v1.R RUN
```

仅将public/同步到本结果目录，再在有pandas/matplotlib/openpyxl环境运行`python brca_malignant_wgcna_present_v1.py RESULT_DIRECTORY`。网络无需任何LYPLA1专用调参。

首次准备因Pal缺少8个Wu同名基因被严格检查终止；随后明确缺项处理并按原种子和参数重跑，未因结果改变构网参数。错误日志留服务器。当前版本提供完整缺项表。

文件入口：BRCA_LYPLA1_WGCNA_figures.pdf、BRCA_LYPLA1_WGCNA_results.xlsx；全部数值TSV保留原WGCNA字段。共有统计字段的对应：cancer=BRCA、stage=06_EXTERNAL、run=20260923T043651Z_malignant_wgcna_v1；kME是bicor、GO效应为overlap并配hypergeometric P/BH q、preservation是Z统计量，无共同的“重要性P值”。不适用的区间/P/q不补造。

方法参考：[WGCNA](https://doi.org/10.1186/1471-2105-9-559)、[hdWGCNA的metacell方法说明](https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html)。后者只作方法参考，本轮实现有上述明确差异。
