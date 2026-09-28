# PRAD LYPLA1 表达 UMAP

## 本轮问题
展示 LYPLA1 在 PRAD 单细胞中的表达位置，并与细胞类型及组织来源对照。

## 输入与范围
复用 PRAD24 CELLxGENE 版本 68b23fda-7191-46a5-8870-819feca3e66e 的 X_umap，全部 68,322 个细胞：癌组织 39,755、癌旁 28,567。不重算降维、不重新聚类、不抽样。原始数据和逐细胞绘图数值留 server165。

## 实际结果
- LYPLA1_UMAP_cancer.png/.pdf：癌组织，左细胞类型，右 LYPLA1 表达。
- LYPLA1_UMAP_all_cells.png/.pdf：全部细胞的注释与表达对照。
- LYPLA1_UMAP_tissue_comparison.png/.pdf：癌旁与癌组织使用相同坐标范围和色标。
- plot_group_summary.tsv：绘图细胞的 pooled 描述汇总，不替代此前患者等权统计。

## 新手解释
颜色越深表示单细胞 log1p(CP10K) 表达越高，灰色是零表达。色标上限为全部阳性细胞表达的第99百分位；超过上限的细胞仍保留并使用最深颜色，箭头表示截顶。较高表达细胞绘制在上层，因此密集群的视觉颜色不等于群体平均表达。

注释取作者 celltype_major_v2；上皮按作者 malignant_anno_merged 拆分，合并显示标签标注 Author-derived，不是新增恶性预测。编号与右侧图例对应，虚线仅为主细胞云的显示轮廓，不是新的亚群边界。

## 限制与反证
UMAP不提供恶性/正常上皮差异检验P值；簇间距离、形状不用于推断功能或轨迹。原始UMAP参数未重新验证。相同色标便于展示，但细胞密度、组成及供者差异仍影响视觉判断。

## 当前决定
本次为表达可视化 DONE，差异检验 NOT_RUN；不根据图像给出显著性或因果结论。

## 下一步
如需判断恶性与正常上皮是否显著不同，另做明确统计单位的表达比较。

## 复现
在 server165 新建独占运行目录（含 private 和 public/06_EXTERNAL），放置 code/prad_lypla1_umap_v1.py；使用 /home/xuzx/miniconda3/bin/python 执行脚本 --out <新运行目录> --commit <代码提交号>。源文件固定来自同级 20260922T142000Z_discovery_v1/source/PRAD24_cellxgene.h5ad。参数、软件版本、文件哈希和细胞数量见 analysis_spec.json、source_manifest.tsv、validation.json。