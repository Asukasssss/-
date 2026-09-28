# COAD 全35基因细胞来源扩展 v2

## 本轮问题

按用户“都计算”的要求，将已经完成的HDC/GSTA4来源核查扩展至本批全部35个功能核查基因；范围不是完整458基因计划池。

## 输入与范围

Lee GSE132465：63,689细胞、23患者、33标本，33,694基因；Pelka作者同步导出：371,223细胞、62患者，43,078基因。使用与v1完全相同的完整矩阵和元数据，SHA256逐项一致。按研究、组织、富集、技术、原始注释及患者分组。锁定方案提交bece086fbc0964537071eaa7c5bd93f4cb649e89先于新表达计算。

原35基因名单逐一保留；新增提取33个基因，HDC/GSTA4汇总和全基因UMI归一化分母精确复用v1。每患者至少20细胞、至少5名达标患者为描述性覆盖规则，同时保留全部覆盖患者；少于3患者不公开表达值。

## 实际结果

两研究各34/35个基因可用；GSTT2均无独立条目，保留不可评估而非填0。ADSS2使用经HGNC确认的旧名ADSS；没有用GSTT2B替代GSTT2。

共85,820行：Lee5,180行，Pelka80,640行；每基因2,452行，按字母顺序分35份TSV。每行保留实际患者数、细胞数、检出率与患者伪合并CPM分布、不可评估原因。HDC/GSTA4旧字段逐字符串一致；本轮没有新P/q，也没有重算患者协变量。

完整35基因的紧凑展示见[全35基因中文总览](ALL35_OVERVIEW_CN.md)，统计数值见[tumor_lineage_overview.tsv](tumor_lineage_overview.tsv)，各基因全部分层见by_gene/。该总览是同一批数值的展示，不是新增独立证据。

在肿瘤、未CD45富集、覆盖达标的谱系中，AQP9、PLA2G15、SLC15A3、SLC29A3、SLC7A8在两研究/技术分层的最高患者中位CPM均位于髓系；BCAT2、CANT1、GSTO2、UCKL1均位于上皮；B4GALT2、BPGM、GNPAT、GSTA4、GSTT2B、ICMT、NNMT、SLC29A1均位于基质。这里的“一致”仅指描述性最高值所在谱系一致，不是显著性或特异性结论。

KMT5A、SLC6A6、PRMT7等的最高值谱系随研究/技术变化，不能强行指定单一来源。ENTPD3、PGAM2、SLC38A3等即使有最高值，其检出率仍很低，需同时查看实际数值。UPP2在达标谱系中的患者中位数均为0，不等于所有患者与细胞均无表达。HDC的Pelka肥大细胞定位与Lee覆盖不足的限制均精确保留。

## 新手解释

这里回答“哪些被采样的细胞类群能检测到这些基因”，并按患者等权汇总，避免细胞多的患者支配结果。CPM是表达尺度，不是酶活、代谢通量或细胞对组织代谢物的贡献。最高中位数只是覆盖足够的类群之间的描述；不是该基因只能来自该类群。中位数为0仍可能在部分细胞或患者检测到。

## 限制/反证

CD45富集、混合富集与未富集分别分析；不计算真实组织细胞比例。Lee经Ficoll处理，Pelka温和解离对基质回收有限。Lee肿瘤肥大细胞覆盖不足，不能把HDC的两研究结果称为稳健重复验证。GSTT2缺条目是测量限制，不是阴性。PRMT7为基因级表达，仍不能提供论文变体比例或PSI。原始上皮注释不等于已确定的恶性细胞。

本批无配套单细胞代谢物，不能验证CAMP代谢物—RNA关系。原M1/M2全检验族均无FDR支持的结论保持。技术、纯度和外部代谢关系验证仍未完成。

## 当前决定

35基因来源表达批次已完成，外部验证与最终整合阶段仍PARTIAL。保留完整974条、原39条及9优先/21保留/5暂挂安排；不按来源表达强弱重新排名靶点。不将未富集总览之外的类群排除，完整富集与亚型结果均可追溯。

## 下一步

结合35基因完整分层与已有功能条件，分别选择适用细胞模型。低检出基因先评估覆盖和独立测量，不能靠来源图补成阳性。HDC仍需适用肥大细胞/空间核查；GSTA4需区分上皮与基质条件；PRMT7变体问题按已有探针边界处理。若开展细胞组成、空间或代谢关系检验，另锁新版本。

## 复现命令

服务器新建独占运行目录并创建.running，复制代码与analysis_spec.json；复用的v1文件须保持原路径及哈希。个体/细胞级文件留server165，导出仅限public_Lee与public_Pelka。

```bash
python3 run_cell_expression35_v2.py --data-dir "$DATA_DIR" --out "$NEW_RUN_DIR" --spec analysis_spec.json --lock-commit bece086fbc0964537071eaa7c5bd93f4cb649e89 --study Lee
python3 run_cell_expression35_v2.py --data-dir "$DATA_DIR" --out "$NEW_RUN_DIR" --spec analysis_spec.json --lock-commit bece086fbc0964537071eaa7c5bd93f4cb649e89 --study Pelka
python code/coad/export_cell_expression35_v2.py --input PUBLIC_AGGREGATES
python code/coad/test_cell_expression35_v2.py
```

源码：code/coad/run_cell_expression35_v2.py、cell_origin_v1.py、export_cell_expression35_v2.py。输入哈希、软件版本和执行校验见source_manifest.tsv与validation.json。两项合成回归检验、实际85,820行唯一键/分位数/公开阈值/精确复用核验已通过。
