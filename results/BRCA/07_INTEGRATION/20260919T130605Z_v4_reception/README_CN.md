# BRCA 117基因 v4 回接结果

## 本轮问题

把文献证据追加到完整历史比较表，形成可追溯的新版本。

## 输入与范围

117个基因、174条关系；基线与交付包详情见 [接收报告](../../05_FUNCTION/20260919T130605Z_v4_reception/README_CN.md)。

## 实际结果

完整表保留原54列，追加47列文献字段及8列接收备注，共109列。117行原顺序和6,318个旧单元格均一致。没有新患者检验或DepMap分数。

- [新手先看：117基因简明表](candidate_reading_view_117.tsv)
- [完整109列比较表](candidate_comparison_117_receiver_v4.tsv)
- [文献编号检查122行](bibliographic_id_check_122.tsv)
- [实际校验结果](validation.json)

## 新手解释

功能文献分类说明哪些对象已有适用研究；工作优先级说明下一项工作安排。两者均不能当作药物靶点认证。全部117基因都保留。

## 限制／反证

接收端核查了编号与选定来源，没有逐篇核实全部实验。5条更正/编辑说明中3条内容影响仍待核定。UPP1需保留全身敲除与细胞特异机制的区别。具体出处和决定见接收报告。

## 当前决定

DONE仅指本次完整表回接；功能取证深度仍PARTIAL。读取最新决定时结合 `review_v4_` 与 `receiver_v4_` 字段，历史列原样保留。上传到 `analysis/brca-functional-review-20260919`，不声称main已合并。

## 下一步

保留24/57/36工作安排，依适用模型与可用数据推进；不把52个文献类别对象当作52个已证实靶点，不把旧10个相关显著基因当作唯一入口。

## 复现命令

```bash
python code/brca_receive_functional_v4.py --repo-root . --output-root /absolute/path/to/new_receipt_check
```

输入参数和来源清单共用05_FUNCTION接收目录，避免重复存成独立证据。
