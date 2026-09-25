# GBM全部候选的两研究细胞表达来源

## 本轮问题
142候选基因主要在哪些细胞类别表达？

## 输入与范围
Darmanis2017原始3589细胞/4供者标签；GBmap内仅Neftel2019子集19030细胞/23供者标签，排除图谱中Darmanis重复数据。沿用作者/图谱注释，不重聚类。

## 实际结果
{
  "status": "DONE",
  "candidate_genes": 142,
  "eligible_profile_rows": 968,
  "evaluable_rank_rows": 233,
  "same_top_all": 45,
  "shared_top_same": 96,
  "no_P_q_tests": true,
  "full_gene_library_normalization": true,
  "donor_equal_weighting": true,
  "all_genes_retained": true
}

## 新手解释
先在每供者每类别汇总，再等权平均。至少20细胞、3来源供者才展示；灰色/NA表示不可评估。1000次整供者重采样的最高类别保持率仅为描述。

## 限制/反证
两研究注释粒度、组织取样与归一化背景不同；全部原取样区域一起描述，不称纯肿瘤核心。跨研究供者别名未独立核验，不称代谢轴独立验证。最高表达不等于唯一作用细胞；没有酶活、通量或机制推断。

## 当前决定
保留全基因及缺测记录，按共同粗类别和各研究全部类别分别比较。

## 下一步
回接关系级/基因级比较表，生成全候选点图与热图。

## 复现
先运行gbm_extract_sc_v1.py，再运行gbm_sc_source_v1.py --out AUDITED_RUN --code-commit COMMIT。
