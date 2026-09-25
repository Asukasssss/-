# PRAD 髓系细胞 inferCNV

## 问题

检查作者注释的髓系细胞是否出现 RNA 推断的染色体偏移，并观察癌与癌旁是否不同。此次运行不把髓系细胞自动改标为恶性。

## 输入

- 同一 PRAD24 单细胞研究，CELLxGENE 数据集版本 `68b23fda-7191-46a5-8870-819feca3e66e`，10x 3′ v3 原始整数计数。作者 major/minor 注释不变。
- 24 位患者中 22 位满足至少 10 个髓系细胞及足够 T 参考细胞的条件，纳入 2,720/2,743 个髓系细胞：癌 1,548、癌旁 1,172。排除的 23 个髓系细胞均来自癌组织。
- 对照包括 B 细胞 622 个、作者标注恶性上皮 1,095 个，以及两组互不重叠的 T 细胞 3,362/3,368 个；每次运行包含相同的选定细胞。总计 11,167 个细胞。
- 每位患者独立运行，参考 A/B 交换，共 44 次。T 参考排除作者 NK、NK-like 和 IFN 类，按组织分层抽样；不使用髓系细胞作参考。
- GENCODE v36 / GRCh38 常染色体基因顺序，31,198 个唯一 Ensembl ID 精确匹配；这是选用的排序参考，不是对原始比对版本的新验证。

## 结果

44/44 次完成，所有纳入细胞均出现在两套参考结果中。各次表达过滤后保留 2,910–8,583 个基因。所有输入细胞常染色体计数大于 100；髓系最低 471、中位数 6,905。

**髓系总体没有表现为癌侧一致增强；当前结果不足以把髓系细胞判为恶性。**

下表为主参考下，每位患者先计算细胞平均偏离分数，再取患者间中位数。偏离分数为 `mean(abs(inferCNV ratio - 1))`，是描述性量，不是 DNA 拷贝数、阳性率或恶性概率。

| 群体 | 癌侧分数 | 癌旁分数 | 同患者癌侧更高/可配对人数 |
|---|---:|---:|---:|
| 髓系总体 | 0.016336 | 0.015719 | 8/15 |
| 巨噬细胞 | 0.018095 | 0.017091 | 10/14 |
| 单核细胞 | 0.015932 | 0.015730 | 7/13 |
| cDC2 | 0.013656 | 0.013846 | 5/15 |

癌/癌旁列分别涵盖所有可用患者，人数可能不同；最后一列才是同患者比较。髓系总体在第二参考下仍为 8/15 对升高；巨噬细胞仍为 10/14、cDC2 仍为 5/15。单核细胞变为 6/13，方向不稳。这里没有进行新的 P/q 检验，不把描述性差异称为显著。

主参考对照分数：留出 T 细胞 0.011852，B 细胞 0.014438，作者恶性上皮 0.023165。上皮对照分数高于髓系总体，但这些群体之间的差异不能单独作为恶性判据。

同一髓系细胞在两参考下的 10 Mb 图谱相关性较高，各亚群×组织的相关系数中位数约 0.922–0.960。这只说明对这两组随机 T 参考较稳定，不是独立验证，也不能排除共同的细胞类型表达偏差。

## 解释

总图中髓系的 chr6:30–40 Mb 红带在癌和癌旁都出现：总体平均推断 log2 ratio 分别为 0.179780、0.182059。坐标表显示该区包含 HLA-DRA、HLA-DRB1、HLA-DPA1、HLA-DPB1；需要警惕免疫表达差异造成的信号，不能把这条共同红带直接当作肿瘤特异 DNA 扩增。本次没有证明该信号的具体机制。

巨噬细胞在较多配对患者中轻度升高，可保留为探索观察；单核细胞基本接近，cDC2 多数反而降低。癌减癌旁图中，髓系总体及主要亚群的差异较弱，不能支持“癌内髓系普遍具有更强 CNV”的表述。

## 限制

- inferCNV 从 RNA 表达推断信号，受谱系、细胞状态和参考选择影响；没有 DNA 层验证。作者恶性上皮是比较对象，不是独立 DNA 验证。
- HMM=False：没有离散扩增/缺失片段调用，没有新增恶性标签，没有预设或套用恶性阈值。原髓系 `malignant_anno_merged=not_applicable` 保持不变。
- 只分析常染色体；按表达过滤后各患者基因集合不同。10 Mb 汇总窗口至少含 10 个保留基因，窗口可用人数不同，具体见表；灰色表示不可用，不是阴性。
- 配对描述要求同患者两侧存在该亚群，没有额外设置每侧亚群细胞数门槛。pDC 只有 1 对可作同亚群配对比较；mDC_regs、cDC1、肥大细胞等部分组细胞少。醒目的稀少亚群色带不能替代足够的患者支持。
- 主图是患者等权的亚群汇总，可能掩盖稀有异常细胞；逐细胞分数、图谱及最终 R 对象留在 server165，不能据汇总图排除所有稀有异常。
- 同研究两套参考不是独立队列。患者身份沿用作者元数据；基因型匹配为 NOT_RUN（无基因型数据）。

## 决策与下一步

不根据当前结果创建“恶性髓系”标签，也不据此重分 SLC6A6 高低组。若进一步追查，优先针对巨噬细胞的具体异常候选核查双细胞、上皮标记污染和连续大片段信号，并以独立 DNA 或其他正交证据验证；这部分尚未执行。

## 文件与复现

- `cohort_infercnv_heatmap.png/.pdf`：两套参考的患者等权热图。ALL_MYELOID 包含下方亚群，行间不可相加。显示范围 ±0.25 仅用于着色。
- `paired_tumor_adjacent_profiles.png`：同患者、同亚群癌减癌旁图；显示范围 ±0.1 仅用于着色。
- `inferred_deviation_scores.png`、`group_score_summary.tsv`：患者均值的中位数和四分位区间。
- `paired_tumor_adjacent_score_summary.tsv`、`paired_tumor_adjacent_bin_summary.tsv`：配对方向及窗口差异。
- `cohort_bin_profiles.tsv`、`group_profile_burden_summary.tsv`：窗口均值、覆盖人数及描述性强度。
- `reference_sensitivity_summary.tsv`、`myeloid_cell_reference_sensitivity_summary.tsv`：同群体/同细胞参考稳定性。T reference/holdout 在两次运行中交换成员，不计算这两个非匹配群体间的稳定性相关。
- `analysis_spec.json`、`source_manifest.tsv`、`code_manifest.tsv`、`conda_explicit.txt`、`R_sessionInfo.txt`、`infercnv_run_defaults.txt`、输入/完成验证文件：来源、参数与版本。

inferCNV 1.26.0，R 4.5.3；cutoff=0.1、min_cells_per_gene=3、window_length=101、cluster_by_groups=TRUE、analysis_mode=samples、denoise=TRUE、sd_amplifier=1.5、HMM=FALSE。未额外按文库大小删除选定细胞；输入检查确认全部超过包默认 100 计数阈值。

服务器运行根目录：
`/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/PRAD/B/20260925T141451Z_myeloid_infercnv_v1`

先在新的独占目录创建 `source/`、`private/`、`public/06_EXTERNAL/`，复制本次五个 `code/prad_*myeloid_infercnv_v1.*` 脚本。用 `prad_prepare_myeloid_infercnv_v1.py --out <新目录> --source-run <20260922T142000Z_discovery_v1目录> --code-commit <代码提交>` 准备输入。在 server165 本地盘用 `conda_explicit.txt` 建独立环境，将路径写入 `source/env_path.txt`。随后依次运行 `prad_qc_myeloid_infercnv_v1.py --root <新目录>`、`prad_schedule_myeloid_infercnv_v1.py --root <新目录> --pilot`、去掉 `--pilot` 的批次命令，以及 `prad_summarize_myeloid_infercnv_v1.py --root <新目录>`。保留私有中间结果，仅交付 `public/06_EXTERNAL/`。

Python 脚本使用 `/home/xuzx/miniconda3/bin/python`，相关版本见 `python_runtime.json`；R 由调度脚本调用独立环境中的 Rscript。

初次共享盘环境安装出现文件复制权限错误，改用服务器本地盘后成功；首例 Matrix Market 稀疏格式兼容问题通过转换成 CsparseMatrix 修复。失败日志保留，计数值和分析阈值未改。正式 44 次运行全部完成；绘图调整没有重跑 inferCNV。

方法参考：[官方 inferCNV 运行文档](https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV)。原始数据、逐细胞/患者测量与运行对象均仅保存在 server165。
