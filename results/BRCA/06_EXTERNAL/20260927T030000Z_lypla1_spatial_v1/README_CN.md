# LYPLA1：Wu2021 BRCA空间分析第一轮

## 本轮问题
检验LYPLA1是否在作者病理标注的含癌区域富集，并评估样本间重复性。主比较为含癌区域与Stroma；其他区域比较为探索性。

## 输入与范围
Wu2021 Zenodo4739739过滤计数矩阵和作者元数据；复用同研究CELLxGENE H5AD的坐标和H&E图像。六份样本对应六个不同的作者patientid（4 TNBC、2 ER）。矩阵与元数据MD5通过，所有条形码、总UMI、检出基因数及H5AD病理标签一致。15611个spot，排除Artefact、Uncertain及缺失标签后15324个用于合适区域的比较。源数据和逐患者/spot数值留server165。

## 实际结果
LYPLA1在9236/15611个spot中检出（59.16%；各样本38.69%–88.82%）。这是覆盖统计，分母包括所有过滤矩阵spot。

|比较（每区域至少20个spot）|患者数|含癌区域较高|患者层面双侧P|解释|
|---|---:|---:|---:|---|
|含癌区域 vs 间质|6|6/6|0.03125|支持区域富集，差异较小|
|含癌区域 vs 混合正常区|2|2/2|不作正式检验|数量不足，且非纯正常上皮|
|含癌区域 vs 明确正常导管|0|不可评估|NA|仅一张切片有3个Normal duct spot|
|含癌区域 vs DCIS|1|0/1|不作正式检验|不能判断浸润过程变化|
|含癌区域 vs 淋巴区域|4|4/4|0.125|同向但未达P<0.05|

主比较平均差为+0.1614 log1p(CP10k)，不是倍数变化。P采用患者差值方向的精确符号检验，非合并spot检验。补充Wilcoxon P亦为0.03125。按用户要求不以FDR作显著性门槛。

## 新手解释
这六份乳腺癌样本中，含癌位置平均比间质位置更容易出现较高LYPLA1信号，说明既有单细胞发现具有一定组织空间支持；尚不能说是纯癌上皮对纯正常上皮的验证。

## 限制与反证
- 改用区域总计数/总UMI的pseudobulk归一化，主比较为5/6同向、1/6反向，说明强度及一致性受汇总方式影响。
- 在各样本内用总UMI作协变量后，6/6仍为正，平均调整后差异缩小至+0.0683。该回归仅作效应方向敏感性，不使用spot回归P值。
- 最低区域spot数改为10或30，主比较仍纳入6例且6例同向；这些是同一数据的敏感性，不是独立验证。
- spot含多种细胞；未完成肿瘤细胞比例调整，不能证明单个癌细胞内上调或LYPLA1特异属于癌细胞。
- 两例混合正常区和3个正常导管spot不足以支持正常上皮正式验证。
- 未建立侵袭边界或细胞邻接模型；没有检验促侵袭、免疫招募、酶活或代谢通量。
- 同研究的单细胞与空间可能患者重叠，不增加为独立队列。
- bootstrap区间仅为小样本描述；显著性以预先指定的精确检验P为准。

## 当前决定
第一轮表达覆盖、作者病理区域比较和六张空间图已完成。结果是有限但一致的空间支持，不能包装成强机制结论。

## 下一步
优先核对作者提供的细胞丰度估计，评估是否可控制癌细胞含量；随后在Andersson队列做同样比较。以上未在本轮执行。

## 复现命令
环境：server165 /home/xuzx/miniconda3/envs/ov/bin/python。主计算代码提交2abec1e。

```bash
python brca_lypla1_spatial_v1.py --root /public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/BRCA_Wu2021_Visium_Zenodo4739739 --out /public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A/NEW_RUN_ID --commit CODE_COMMIT
python brca_lypla1_spatial_audit_v1.py /public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A/NEW_RUN_ID
```
先建立新运行目录及独占.running锁。坐标文件复用路径及可下载版本见source_manifest.tsv和spatial_source_citations.json。患者和spot表、六张空间图、详细输入哈希保存在本轮服务器目录；本仓库仅保留队列汇总和代码。
