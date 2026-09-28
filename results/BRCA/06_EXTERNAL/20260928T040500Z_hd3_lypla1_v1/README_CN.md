# LYPLA1：乳腺癌 Visium HD 3′ 空间探索

## 本轮问题
寻找 HD 空间转录组并实际分析 LYPLA1 的空间表达；优先保留癌区/正常上皮比较所需的证据边界。

## 输入与范围
采用 [PacBio 公开数据](https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/DATA-RevioSPRQ-Kinnex-VisiumHD-humanBreast/)，Visium HD 3′ 建库联合 Kinnex 长读长测序，供者为一位70岁女性，来源说明为浸润性导管癌、新鲜冰冻组织。技术来源：[厂商应用说明](https://www.pacb.com/wp-content/uploads/Application-note-Kinnex-single-cell-RNA-kit-for-10x-Visium-HD-3-Spatial-Gene-Expression.pdf)。分析使用基因级 Kinnex 处理矩阵，不将报告中的短读长流程统计冒充本矩阵统计。

输入为22,854个特征、3,672,073个2 μm空间条码、18,721,468个非零项。原始 LYPLA1 特征计数共11,408；另外 LYPLA1+TCEA1 联合特征的1个计数未并入目标基因。

另核查10x常用探针型 HD 乳腺 DCIS 数据：LYPLA1 在参考表被标为 included=FALSE，官方 feature_slice 的 target_sets 也未保留它。即使原始层有读数，也未用其作可靠表达证据，见 probe_coverage_audit.json。此结论针对核查的v2探针及数据，不推及所有HD技术。

## 实际结果
|网格|合格网格|LYPLA1 检出网格|检出率|
|---|---:|---:|---:|
|8 μm|155,877|9,820|6.3%|
|16 μm|46,717|8,486|18.2%|
|32 μm|15,322|5,659|36.9%|

16 μm为预设主展示尺度，8与32 μm为敏感性展示。网格增大导致检出率上升，不表示基因被进一步上调。

LYPLA1 与上皮标记均值的描述性 Spearman 相关：8 μm为0.126，16 μm为0.257，32 μm为0.363；控制总计数的偏秩相关分别降为0.019、0.042、0.057。完整结果见 marker_associations.tsv。

空间覆盖及可视化完成；正式“癌上皮高于正常上皮”仍为 NOT_EVALUABLE，没有基于网格数量计算显著性 P 值。

## 新手解释
这次确实在HD数据中测到了LYPLA1。较粗尺度能让分布更容易看清，但其与上皮标记的同步变化很大程度随测序深度变化，校正后关联较弱。不能把这些图解释为癌细胞特异富集或癌上皮高于正常上皮。

## 限制/反证
仅一份标本；公开下载包和web_summary未提供明确病理癌区/正常导管标签。EPCAM/KRT等标记不能区分恶性与正常上皮。HD的2 μm条码及8/16/32 μm网格都不是已分割的单细胞；稀疏零值不能当作生物学不表达。上皮关联经深度校正后减弱是必须保留的限制。

## 当前决定
保留为单标本HD空间探索，不升级为正常上皮对照验证。没有继续强制用被过滤探针作阳性图。

## 下一步
若要回答癌上皮高于正常上皮，需获得这一标本可追溯的病理区域标注及足够正常导管，并在更多独立供者重复。当前图可用于查看位置和规划后续标注。

## 复现命令
本地 `python code/brca_hd_fetch_v1.py` 将源文件经RAM传到server165，不保存本地源矩阵。服务器执行 `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/xuzx/miniconda3/envs/ov/bin/python brca_hd3_analyze_v1.py RUN_DIR`。随后用 `brca_hd3_zoom_v1.py RUN_DIR` 从保存的16 μm结果绘图，不重跑统计。分析代码提交 bde0827。源文件 SHA256全部核对一致。

图像对齐采用来源web_summary的变换矩阵；本矩阵8 μm总计数与提供方UMI图的秩相关为0.961，组织掩膜427,917个点与原报告完全一致。坐标变换后的QC点均在图像范围内，另完成图像叠加目视检查。
