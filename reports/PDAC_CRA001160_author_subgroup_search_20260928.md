# CRA001160 原作者七亚群标签检索

日期：2026-09-28。结论：找到原作者七亚群的标记基因表，但本轮检查的公开来源中尚未找到逐细胞的七亚群归属。不能据此断言作者从未公开；也不能把自行聚类当作作者标签。

## 已实际检查

1. [原论文](https://www.nature.com/articles/s41422-019-0195-y)及 PMC 正文、附件索引：数据可用性仅指向 CRA001160/PRJCA001063；未找到独立代码仓库链接。附件为6份补充图PDF及7份XLS表。
2. [GSA公开目录](https://download.cncb.ac.cn/gsa/CRA001160/other_files/)：PDAC.tar.gz、all_celltype.txt、count-matrix.txt。重新列出压缩包140个成员，是35个标本目录及各自genes.tsv、barcodes.tsv、matrix.mtx，没有Seurat/RDS或细分注释对象。
3. all_celltype.txt：两个字段cell.name、cluster，10个大类；导管细胞只有type1/type2。
4. Table S1–S7：均已在服务器下载、读取全部工作表与字段。S1临床资料；S2逐细胞QC；S3大类标记；S4拟时序基因；S5七亚群标记；S6 TCGA使用的标记；S7 TCGA组间差异。没有逐细胞亚群归属表。
5. Besca/Zenodo 3969339：实际读取obs注释；原始Cell_type仍为10大类，衍生leiden为54群，celltype3为17类。其更细标签主要细分免疫等类别，不能当作原作者type2七亚群。
6. GitHub检索CRA001160、all_celltype.txt、41422_2019_195及相关亚群词。检索到的后续研究使用原大类注释或自行重聚类。其中[PDAC_MSLN代码](https://github.com/Oliver-Liang-1999/PDAC_MSLN/blob/691602bba7c524cf2abcbce5f61dbb3eadf520d1/Human%20scRNASeq%20analysis_v2_revision.Rmd)虽引用S5，但明确为自己重新划分的4群，并非原7群的逐细胞标签。

## 原作者 Table S5

[直接下载](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41422-019-0195-y/MediaObjects/41422_2019_195_MOESM11_ESM.xls)

文件：41422_2019_195_MOESM11_ESM.xls；321,536字节；SHA-256：5a708b67b13c143984f606f0d79d87db0e928ea1bdd4c581f8547f75b44026ff。

包含2,282条基因—亚群标记记录、1,910个不同基因，字段为Gene、p value、avg_logFC、pct.1、pct.2、Subpopulation。

|作者亚群|标记记录数（不是细胞数）|
|---|---:|
|1|349|
|2|468|
|3|60|
|4|234|
|5|186|
|6|625|
|7|360|

LYPLA1未出现在此标记清单中。这不表示LYPLA1不表达，也不能据此判断七群间没有差异。该表没有条码或每个细胞表达值，不能计算要求的LYPLA1七群比较。

S2含57,530个细胞的QC记录，字段为Cell ID、Number of genes、Number of UMIs、Percentage of mitochondria gene；没有Subpopulation列。S2及临床数据均留服务器，未回传本机。

## 下一步所需材料

精确复现作者七群需要：
- Fig.3e / Fig.S3 对应的 type2 细胞条码 → subgroup1–7 映射；
- 或包含该标签的作者 Seurat/RDS/其他分析对象；
- 若还有type1两亚群映射，可同时提供；附条码格式与样本前缀说明。

论文公开通讯联系人包括 Yun-Gui Yang（ygyang@big.ac.cn）、Wenming Wu（doctorwuu@126.com）和 Yupei Zhao（zhao8028@263.net）。本轮未向任何人发送消息。

可复制的索取内容：我们希望使用CRA001160原始细胞归属比较LYPLA1表达。请问能否分享论文Fig.3e中type2导管细胞七亚群对应的cell barcode/sample ID/subgroup映射表，或保留该注释的Seurat对象？如有type1两个亚群的映射，也希望一并获得。我们不需要重新索取原始测序数据。

## 分析状态

已停止此前自行聚类进程并保留其输入与中间文件；没有将其标为作者七群分析完成。本轮没有产生作者七群的LYPLA1差异数值。原作者标签未取得前，任何依据S5标记进行的映射只能另列为推断性注释，不能称为原作者逐细胞标签。
