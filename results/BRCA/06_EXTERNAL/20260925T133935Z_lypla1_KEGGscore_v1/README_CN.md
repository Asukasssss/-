# LYPLA1／LYPLA2：KEGG甘油磷脂代谢评分补充分析

## 本轮问题

按用户要求更换为一个预先固定、范围更广的基因集，并加快执行。新增KEGG hsa00564甘油磷脂代谢，完整保留上一套Reactome结果；不以显著性决定最终展示哪一套。按名义P保留探索线索，q存档。

## 输入与范围

来源为[KEGG官方hsa00564](https://www.kegg.jp/entry/hsa00564)，API与获取时间／数据库说明见`kegg_provenance.json`。105条GENE记录中，一条hsa:128966744没有符号，仅有KO/EC，明确标为UNRESOLVED，不当作基因名称或推测对应。

按上一轮相同覆盖规则，固定**92个共同可评估成员**，排除LYPLA1、LYPLA2以及身份不明／缺测／低检出成员。完整105条仍保存在`kegg_members.tsv`与`frozen_score_members.tsv`。未按相关性挑基因。

- 与旧24基因集合重叠20个；新增独有72个，旧集合独有4个；Jaccard=20/96=0.2083。这不是简单更换集合名称。
- Wu仍为20位供者、24,489个恶性上皮细胞；Pal仍为33位供者、130,078个恶性上皮细胞。Pal沿用Chen2026重注释，不称原Pal标签。
- 每供者原始counts汇总后log1p(CP10k)；三个队列分别按基因z标准化，再将同一92成员等权平均。CAMP保持作者RNA尺度后标准化。不同队列分数绝对值不可直接比较。
- CAMP主分析60例，作者可用值55例，组织类型冲突记录继续排除。
- 完全沿用8项检验：两套单细胞×LYPLA1／LYPLA2=4项，CAMP主分析和可用值=2项，同一ER＋组成PC12调整=2项。99,999次置换、4,000次bootstrap、独立供者／病例为单位，参数未为结果优化。
- 组成代理按同一规则排除156历史候选和本轮完整通路成员，实际仍只移除与标记重叠的KYNU、HAL，9个代理面板；这不是细胞比例或完整纯度估计。
- 仅为避免误读，缺符号记录在任何关联计算前从KO文本改为明确的未解析KEGG编号。选中的92成员和数值没有变化，初始脚本／输入在服务器保留，修正记于`symbol_label_correction.json`。

## 实际结果

|检验|n|ρ／偏秩ρ|逐项95%区间|P|
|---|---:|---:|---|---:|
|Wu：LYPLA1—评分|20|+0.200|−0.277～+0.632|0.39815|
|Wu：LYPLA2—评分|20|+0.305|−0.110～+0.646|0.18898|
|**Pal：LYPLA1—评分**|33|**+0.434**|**+0.076～+0.700**|**0.01141**|
|Pal：LYPLA2—评分|33|+0.215|−0.151～+0.525|0.22567|
|CAMP：评分—LPC主分析|60|−0.159|−0.388～+0.085|0.22372|
|CAMP：评分—LPC可用值|55|−0.145|−0.393～+0.118|0.28850|
|CAMP：主分析＋ER/组成PC|60|−0.054|−0.337～+0.238|0.68924|
|CAMP：可用值＋相同调整|55|+0.006|−0.295～+0.306|0.96634|

仅Pal的LYPLA1—评分达到本轮名义P<0.05。完整新8项在`results.tsv`；旧8项与新8项并列在`both_gene_sets_results16.tsv`，旧效应/P/q全部复用，不重算覆盖。两套各自的BH8留档，不能认为这校正了项目全部历史选题。

## 新手解释

**Pal中LYPLA1与更广泛的甘油磷脂相关表达背景相联系**，可以按用户当前的P标准保留。但Wu尚无明确统计支持，仍不是两个队列稳定验证。

换集合后，Pal的点估计从0.385变为0.434，并不证明新集合更优。两次使用同一批供者，成员存在重叠，也没有正式比较相关系数差异。

**CAMP的新评分没有显示与实测LPC的明确联系**。主分析点估计为弱负向，背景调整后接近零。因此这次同样没有建立“LYPLA1→脂质评分→LPC”的解释链。调整后减弱也不能据此断言全部由组织组成造成。

LYPLA2没有明确名义P支持，但没有直接比较LYPLA1与LYPLA2的关联强度差异，不能把LYPLA2当成已经验证的阴性对照。

## 限制／反证

这是看过上一轮结果之后的补充探索，不是独立验证或替换旧结果。分数是覆盖子集的RNA汇总，既有合成也有分解等不同环节，不代表脂质含量、净代谢活性或反应通量。

细胞汇总使用全部合格恶性上皮供者，不把细胞当独立病例；同时也不能推断每位患者内部LYPLA1高的单个细胞均有此状态。CAMP组织分数不能称为恶性上皮细胞内分数，三个队列未进行同患者拼接。

精确符号匹配不强行使用未核实别名。除两个比较基因外，被排除的11项为PLA2G4B、PLA2G10BP、CHAT、CHKB、未解析hsa:128966744、DGKK、LPCAT4、PLA2G2E、TAFAZZIN、PLA2G12B、JMJD7-PLA2G4B；逐队列原因见`gene_coverage.tsv`。

## 当前决定

保留“Pal中LYPLA1具有甘油磷脂相关转录背景”的线索，同时明确Wu未支持、CAMP评分—LPC未支持。已有LYPLA1—LPC直接关联及配对差值分析保持不变。本轮不足以升级为代谢机制，也不取消原候选。

## 下一步

停止本轮基因集切换，将两套固定结果一起讨论。不继续用更多集合／算法寻找显著性；若另提问题，另行记录假设与范围。

## 文件、核查与复现

- `LYPLA_KEGGscore_results.pdf`：新8项三张图；`two_gene_sets_comparison.pdf/png`：两套全部结果并排。
- `results.tsv`、`both_gene_sets_results16.tsv`：完整新／旧统计；`gene_set_overlap.tsv/json`：成员交集。
- `analysis_spec.json`、`source_manifest.tsv`、`validation.json`、`independent_R_audit.*`：参数与审计。R独立重算三个队列评分和8个相关系数／BH8；没有重复全部置换或bootstrap抽样。
- Excel为本机／服务器便读副本；GitHub提供全部对应TSV、图件、脚本与解释。原始矩阵、逐供者／病例值、身份映射均在server165。
- KEGG原始API文本保留服务器`source/kegg_hsa00564_source.txt`，其SHA256在manifest及provenance中；公共成员表可直接查看。

在新独占运行目录中复制对应配置、固定成员与源文件索引、复用帮助脚本，执行：

```bash
python3 brca_lypla1_KEGGscore_prepare_v1.py "$RUN_DIR"
python3 brca_lypla1_KEGGscore_analyze_v1.py "$RUN_DIR"
Rscript brca_lypla1_PCscore_audit_v1.R "$RUN_DIR"
python3 brca_lypla1_score_comparison_v1.py "$OLD_PUBLIC" "$NEW_PUBLIC"
```

审计脚本和私有`PC_score`列名沿用旧版接口，当前数值实际为92基因KEGG评分，公开通路ID和图名已明确更新。两份新计算脚本基于提交899f9d9的旧脚本，仅适配新成员／版本／名称及无符号记录，不改变统计算法。

server165运行目录：`/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A/20260925T133935Z_lypla1_KEGGscore_v1`。上传分析分支不表示已合并main。
