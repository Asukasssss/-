# PYCR1直接干预数据可用性：GSE220931

## 本轮问题
从相关分析推进到直接PYCR1敲低，检查已有乳腺癌模型公开数据能否检验固定E2F及胶原形成程序。

## 输入与范围
已下载GSE220931处理后矩阵及9份GSM元数据到server165新运行目录。58,735行、9列测量；9列均与作者样本标题精确对应。未下载原始测序/质谱，未导出源矩阵。仅开展有界搜索，不代表穷尽全部公开数据。

## 实际结果
矩阵数值全部非负整数，无缺失，基因ID无重复；但作者GEO处理说明称FPKM，没有明确声明该附件为原始计数。论文方法写hg19，GEO写hg38。
GSM6829547标题shP2_2、genotype却写shP3；GSM6829548标题shP2_3、genotype写shP4。不能自行认定为笔误或额外独立构建。
因此正式程序检验状态NEEDS_REVIEW，新增P/q为0。PYCR1附件数值仅作输入检查，组内中位数为ctrl=1450、shP1=924、shP2=353；单位未确认，不能解释成敲低百分比或倍数，不能当归一化差异。

## 新手解释
找到了直接干预实验，不等于现在就能安全计算。必须知道数字的单位和哪些样本属于同一干预。整数看起来像计数，但不能据此覆盖作者说明；分组标题与基因型不一致也不能靠猜测解决。这是数据解释缺项，不是PYCR1生物学阴性。

## 限制/反证
GSE220931来自MDA-MB-231癌细胞，一个研究内的两种敲低构建不是两个独立队列。样本的独立培养/克隆关系尚未核实。论文已有PYCR1功能发现，我们的再分析不应称为首次发现。
CAF论文的PXD018343包括CAF/NF蛋白及c646干预；PXD024746是同位素标记ECM蛋白资料。不能仅因所属论文研究PYCR1，就将其当成PYCR1敲低转录组。GSE171983是肝癌PYCR1干预，不纳入乳腺癌正式复现。GSE196354是CAF/NF比较，不是PYCR1干预。

## 当前决定
保留全部既有117基因结论和上一批数值。该数据的正式计数差异及两个程序检验暂挂，记录具体缺项，不用其他癌种或其他基因干预填补。

## 下一步
需要能明确说明附件计量单位、shP2_2/3构建归属及独立重复设计的作者补充材料或澄清；之后按两构建分别比较对照，并对固定两程序×两对比统一校正。此次未联系作者、未发送消息，未从原始FASTQ重建矩阵。

## 复现命令
python3 code/brca_pycr1_perturb_inventory_v1.py --root "$SERVER_RUN_ROOT"

源下载URL：https://ftp.ncbi.nlm.nih.gov/geo/series/GSE220nnn/GSE220931/suppl/GSE220931_processed_data.txt.gz
元数据：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220931&targ=self&form=text&view=full；各GSM同样接口。
论文：https://www.nature.com/articles/s41419-023-06200-5
CAF：https://www.nature.com/articles/s42255-022-00582-0
蛋白来源：https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD018343 与 https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD024746
