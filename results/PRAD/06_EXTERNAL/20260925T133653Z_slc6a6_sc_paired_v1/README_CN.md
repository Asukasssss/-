# SLC6A6 单细胞供者配对比较

本轮问题：同类细胞中，SLC6A6在肿瘤组织与癌旁组织是否不同？

输入与范围：沿用作者细胞标签及donor_id，复用前轮供者表达汇总；每位供者两种组织对应类别均至少20细胞，至少8对才做推断。所有同名细胞类别预先纳入；另列恶性上皮与癌旁正常上皮的状态比较，不能当作完全相同细胞状态。匹配仅按明确donor_id，不按顺序推断。

实际结果：
- B-cells: n=2, effect=0.0444529, P=nan, up/down/equal=2/0/0, status=NOT_EVALUABLE
- CAFs: n=11, effect=0.0203514, P=0.174805, up/down/equal=8/3/0, status=DONE
- Cycling: n=0, effect=nan, P=nan, up/down/equal=0/0/0, status=NOT_EVALUABLE
- Endothelial: n=15, effect=-0.000345607, P=0.719727, up/down/equal=9/6/0, status=DONE
- Epithelial_altered_benign: n=12, effect=0.00175715, P=0.339355, up/down/equal=8/4/0, status=DONE
- Epithelial_malignant: n=6, effect=0.0186681, P=nan, up/down/equal=5/1/0, status=NOT_EVALUABLE
- Epithelial_normal: n=4, effect=-0.0493607, P=nan, up/down/equal=0/4/0, status=NOT_EVALUABLE
- Myeloid: n=11, effect=-0.00855176, P=0.637695, up/down/equal=4/7/0, status=DONE
- PNS_glial: n=4, effect=-0.00773459, P=nan, up/down/equal=1/1/2, status=NOT_EVALUABLE
- Plasma_cells: n=0, effect=nan, P=nan, up/down/equal=0/0/0, status=NOT_EVALUABLE
- SMCs: n=15, effect=0.00297943, P=0.635498, up/down/equal=5/8/2, status=DONE
- T-cells: n=15, effect=0.00744155, P=0.454285, up/down/equal=8/7/0, status=DONE
- Malignant_vs_adjacent_normal_epithelium: n=5, effect=-0.0578922, P=nan, up/down/equal=0/5/0, status=NOT_EVALUABLE

新手解释：以供者为统计单位，每位供者先计算同类细胞平均表达，再比较肿瘤减癌旁。P来自双侧精确符号置换的Wilcoxon有符号秩检验；展示按P<0.05，q保留不作为入口。图中区间为平均变化的配对bootstrap区间，与秩检验不是同一统计量。

限制/反证：同一研究的补充分析，不是第二队列验证；低覆盖类别不可评估，不作阴性。表达是平均log1p(CP10K)，不是浓度倍数或转运活性。宽细胞类型内部的亚群构成变化、技术因素仍可影响结果。继承原CELLxGENE固定版本和其GEO链接不一致警告。

当前决定：按真实结果报告，不将组织RNA差异套用到某一细胞类型。

下一步：独立研究和亚群构成检查尚未完成。

复现命令：python prad_slc6a6_sc_paired_v1.py --out NEW_RUN --source-run ORIGINAL_RUN --code-commit PARENT；先建立新运行目录。
