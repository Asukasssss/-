# COAD原CAMP P<0.05结果展示

## 本轮问题

按用户要求展示原CAMP代谢物结果中原P<0.05的项目；不使用本轮33对配对P或代谢物—RNA关联P替代。

## 输入与范围

server165原项目results/tables/cohort_effects.tsv的COAD全159行。原设计为37个肿瘤、39个正常标本，原analysis_design明确Unpaired。沿用作者处理值的项目历史Wilcoxon P，字段wilcoxon_p；不是声称该P直接来自原论文附件。通过feature_name及metabolite_key核对全部159项，Hedges g与wilcoxon_fdr分别逐字符串吻合冻结表hedges_g与effect_fdr。输入哈希见source_manifest.tsv。

## 实际结果

原P<0.05共83项：按原g方向40项肿瘤较高、43项较低。原q<0.05仍为73项；83项中10项仅名义P支持（7较高、3较低）。原数值完整保留，不重新计算、反推P或重新BH。

[全部83项中文结果](ORIGINAL_P_LT_005_CN.md)｜[中文TSV](original_P_lt_005_CN.tsv)。original_coad159.tsv保留全159项原字段与数值，results.tsv按公共统计前缀展示同一批结果及筛选标记。

## 新手解释

原P<0.05表示原肿瘤与正常组比较的名义显著性；q是原多重比较校正结果。原g是标准化效应，不是倍数。这里没有33对患者的逐对人数；前一批33对是不同样本范围和设计的新分析，其94项q显著不替代本批83项原P显著。

## 限制/反证

10项P显著但q不显著保持区分，不称通过FDR。原比较为非配对，虽两组包含同一来源患者的组织，不能用这张旧表声称已进行配对推断；该历史设计的限制不在本轮展示中悄悄修正。处理值填补与化学名称身份限制继承原来源，中文翻译不作为新鉴定。未重跑患者矩阵，未核验原论文逐峰身份；没有新统计或独立验证。

## 当前决定

本批原P视图DONE。原159项、73项原q、39条关系、9/21/5安排均保持。输入/数值一致性和筛选数量已核验。

## 下一步

用户可按原P浏览完整83项；比较其他结果时分别标注原P、原q、新配对P/q和肿瘤内关联P/q。

## 复现命令

在server165使用新的独占运行目录，输入均为已有汇总：

```bash
python3 export_original_camp_p_v1.py --source "$ROOT/results/tables/cohort_effects.tsv" --frozen "$FROZEN_EFFECTS" --template "$STATISTICAL_TEMPLATE" --out "$NEW_RUN_DIR"
```

本机仅渲染汇总：

```bash
python code/coad/render_original_camp_p_v1.py --results results/COAD/01_CAMP/20260921T125510Z_original_camp_p_v1
```
