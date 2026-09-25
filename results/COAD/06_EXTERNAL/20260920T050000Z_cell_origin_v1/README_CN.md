# COAD：细胞来源覆盖与表达前规则锁定

## 本轮问题

在不读取 HDC/GSTA4 表达的情况下，确认 Lee/Pelka 的患者、标本、文库、富集与技术条件，先交覆盖；补取实际 PRMT7 平台探针。协变量不重算。

## 输入与范围

Lee GSE132465 的原作者完整细胞注释；Pelka 作者 colon10x_default 的同步元数据、注释、细胞和文库向量。矩阵及逐细胞/患者映射只在 server165。作者加载代码 readDataRobj.m 和 Figure_1.m 明确这些向量的同步对象关系；不是按患者编号猜配对。

GSE132257 是 2 名患者的 10 个制备/技术样本；GSE144735 是 6 名比利时患者的核心、边缘及正常组织，最终 27,414 细胞。二者本批不分析，不另计独立研究；与主队列患者重叠未在个体层面完全核实。

## 实际结果

- Lee：63,689 细胞，23 名患者、33 标本；肿瘤 47,285、正常 16,404。原注释保留六大类及所有亚类。文库 ID 未提供，记 NA；协议为 10x 3′ v2，未报告 CD45 富集，但使用 Ficoll 与冻存。
- Pelka 作者版本：371,223 细胞、62 名患者、100 个 PatientTypeID、181 个 batchID。每个 batchID 的患者、组织、处理和技术标签一致。PatientTypeID 是作者的患者×组织键，不把它解释为额外独立患者或完整的解剖区域记录。
- Pelka 五种处理分别为 unsorted 231,737；CD45pMACS 122,766；mixUnsortCD45MACS 10,670；LiveMACS 4,942；CD45pCD3nCD19nMACS 1,108。SC3Pv2/v3 分开。
- GEO 370,115 个细胞均匹配作者版本；26 个共同元数据字段逐值一致。作者额外的 1,108 个细胞全部属于 CD45pCD3nCD19nMACS。固定作者版本，不能混接 GEO 矩阵。
- 生成 1,226 行覆盖表，不同注释层级是同一批细胞的不同分辨率，不相加。Lee 肿瘤仅 3 个肥大细胞、2 名患者；Pelka 未分选肿瘤肥大细胞在 v2/v3 各有 692/433 个，达每患者20细胞的患者各11/8名。
- GPL16699 实际平台中 PRMT7 有两个 60nt 探针：A_23_P77437、A_23_P77430，序列见 PRMT7_platform_probes.tsv。取得探针这一步 DONE；论文变体区段/接头匹配仍 NEEDS_REVIEW。平台多转录本注释不等于已证明可区分 V1/V2。
- AQP9/KMT5A 等方向条件表原样接收，不计新增独立文献核验。附件 manifest 全部通过。

## 新手解释

现在知道各类细胞覆盖是否足够，还没有得到 HDC/GSTA4 表达结论。Lee 肿瘤肥大细胞极少，因此即使随后未检出，也不能认定肥大细胞没有信号。Pelka 的未分选和富集资料可以分别回答问题，不能混成组织细胞比例。

## 限制/反证

解离、Ficoll、冻存和分选影响回收；Pelka 基质回收有限。没有粒细胞注释不等于没有该生物来源。上皮不自动改标为恶性。外部单细胞表达不验证 CAMP 代谢物—RNA 轴，不代表酶活或通量。PRMT7 仍缺论文构建/变体序列；现有基因级矩阵不能推算 PSI。

## 当前决定

表达前锁定 config/coad_cell_origin_v1.json：每个患者×原细胞类型至少20细胞、至少5名达标患者，作为粗略描述支持标记，非显著性或功效标准；全部低覆盖组保留。按相同研究/组织/富集/技术汇总，患者等权，报告检出率与全基因 UMI 库大小标准化伪bulk CPM。少于3名患者的表达汇总不公开。无检验、无 P/q；原974/39及9/21/5不变。

## 下一步

在本规则远端提交确认后，另建表达运行目录，完整读取矩阵、校验原始计数与身份、按患者计算 HDC/GSTA4。先完成当前数据可支持的定位，再决定是否需要空间验证。PRMT7继续补论文特异序列，不以数据库isoform编号替代论文V1/V2。

## 复现命令

在 server165 的新独占运行目录复制 code/coad/cell_origin_v1.py 与 audit_cell_origin_sources.py，使用 Python3（numpy/pandas）执行；原运行目录不覆盖。

```bash
python3 cell_origin_v1.py --data-dir "$DATA_DIR" --out "$NEW_RUN_DIR" --mode coverage
python3 audit_cell_origin_sources.py --data-dir "$DATA_DIR" --out "$NEW_RUN_DIR"
```

原始入口：[Lee GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465)、[Pelka GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178341)、[作者数据/代码](https://github.com/matanhofree/crc-immune-hubs)、[Pelka 原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8772395/)、[GPL16699](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16699)。来源文件SHA256在source_manifest.tsv。
