# CAMP → 四癌 LYPLA1：统一重绘汇报 v1

## 本轮问题

把已有结果整理为一套可讲清项目逻辑的 PPT，并从图源运行代码统一重绘 UMAP、空间图和定量图。六癌作为 CAMP 发现背景，LYPLA1 只讨论 BRCA、COAD、PDAC、PRAD；GBM、ccRCC 不进入 LYPLA1 主线。

## 输入与范围

- `source_manifest.tsv` 固定 27 张公开汇总表的 Git 提交、路径、URL、SHA256。不是把各分支整个最新结果重新计算一遍。
- `server_figures/server_render_source_manifest.tsv` 记录服务器原坐标、缓存表达和病理图像来源。细胞、患者、spot 测量表没有复制到本地或 Git。
- 所有 P、效应和原置信区间复用；展示阈值为名义 P<0.05。原 q 保留于汇总表，不将其解释为 FDR 控制后的发现。
- 四癌 UMAP 复用既有坐标与入组细胞；BRCA 使用此前分析者生成的上皮 UMAP，其余使用既有作者/缓存坐标。没有新增降维、聚类、推断恶性身份或批次校正。
- BRCA 由原 raw counts 按既有 log1p(CP10K) 公式提取 LYPLA1 供显示；其余复用冻结表达。BRCA 空间采用原区域归属，PRAD 采用原病理标签；不重新配准或补画癌区。

## 实际结果

交付 22 页 `CAMP_LYPLA1_统一重绘版.pptx` 和实际渲染 PDF；标题、说明和页脚为可编辑文本，图表另提供 PNG / PDF 及重绘代码。

|范围|重绘数量|内容|
|---|---:|---|
|定量图|13|六癌覆盖与精确分子对照、候选 RNA、LYPLA1 关系与敏感性、患者比较、检出分解、状态、三个癌种案例|
|UMAP|4|BRCA 28,844；COAD 31,247（含灰色背景 2,747）；PDAC 61,272；PRAD 5,357|
|空间切片|21|BRCA CTA 14 张 / 3 患者；PRAD Erickson 7 张 / 1 患者|

正文依次展示发现框架、四癌 LYPLA1、细胞与空间定位、COAD UCKL1—尿苷、PDAC/PRAD 牛磺酸—SLC6A6、BRCA ASNS/GLS—谷氨酰胺。

统一规则见 `STYLE_STANDARD.json`：白底，红色升高/恶性，蓝色降低/参照，灰色背景，字体与边框统一；相关图实心表示 P<0.05，空心表示未显著。UMAP 点面积 1 pt²、固定随机顺序；表达色标 0–3 log1p(CP10K)，顶部三角表示超上限，各图截顶比例记录在 `render_validation.json`。轮廓仅为展示外包络，不定义新细胞群。空间图均为 H&E—原病理区域—实测表达，不插值、不平滑表达。

完整图册包括全部 21 张切片，正文 BRCA 按样本名字典序取第一张作布局示例（该张没有免疫区，并非癌—免疫比较代表）；PRAD 按字典序取第一张癌/良性均≥20 spots 的切片。选择不使用表达或 P 值。

## 新手解释

可以讲“从 CAMP 关系候选发现 LYPLA1 在四癌具有统计线索，进一步检查表达来源与空间背景”。需要紧接着说明：BRCA、COAD、PRAD 是 RNA 显著升高；PDAC 是 GPC—LYPLA1 关联显著，RNA P=0.0777。不能简写为同一个关系在四癌重复，也不能写成四癌 RNA 均显著。

## 限制/反证

- 不同癌种平台、样本量和检验设计不同。相同色标便于读图，不能消除批次、测序深度或平台差异。
- 部分单细胞结果仅是细胞等权描述，不能替代患者重复性。PRAD 阳性细胞均值几乎相同，主要差异为检出比例。
- BRCA LPC 配对变化主分析显著，可用值敏感性未显著；PDAC 单细胞另一队列及 GeoMx 敏感性没有一致显著支持。
- BRCA 空间为癌区对免疫区，患者层面 P=0.25，不能当作正常上皮对照；PRAD 5 张可比较切片只有 1 位患者。
- 三个案例为癌种代表关系，未检验癌种交互，不能声称癌种特异性机制。
- 本地没有 Microsoft PowerPoint 渲染环境；已用服务器 LibreOffice 实际导出并检查 22 页，PDF 是固定版面参照。空间图册为派生图，不含逐 spot 数值。Erickson 来源按 CC BY-NC 3.0 标注。

## 当前决定

DONE：38 张图完成重绘，22 页完成组装与实际渲染；原统计保留。27 张源表哈希全部匹配。独立目录交付，不覆盖历史结果。代码与小于 5 MB 的成品/单图提交独立分支，公共代码走 PR，不自动合并。

仓库级 `tools/check_repository.py` 报告 13 个继承自父提交的历史文件问题（Excel 扩展名或超过 5 MB），本次不改动这些历史结果。本次新增文件独立按相同大小/格式规则检查；结果记录在 `repository_validation.json`。本目录 `.gitattributes` 禁止 Git 改写冻结文件行尾，保证跨平台源表 SHA256 可复核。

## 下一步

按汇报时长删减页数即可，优先保留候选发现、四癌 RNA 与关系、患者一致性、空间边界和三个代表案例；不因版面删去反证。逐页提纲见 `逐页讲解.md`。完整空间图册与 HTML 预览由 `ppt_package_v1.py` 写入本地可视化目录，不作为新增独立证据。

## 复现命令

环境版本见 `software_versions.json`。本地需要 numpy、pandas、matplotlib、python-pptx、Pillow、PyMuPDF；服务器还需 h5py、scipy。字体 Noto Sans CJK TTC 放入忽略的 `runtime/ppt_fonts/NotoSansCJK-Regular.ttc`；可从 server165 `/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc` 获取。PPT 使用 Noto Sans CJK SC；缺字体时优先阅读 PDF。

```text
# 源表已经随本版冻结；仅重新提取源表时才需要相关 Git 历史对象
python code/ppt_prepare_sources_v1.py
# 在 server165 专属运行目录执行，输出 public/；不在本地加载原始矩阵
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/xuzx/miniconda3/envs/ov/bin/python ppt_server_redraw_v1.py
# 仅把 public/ 派生图及清单复制为本版 server_figures/
python code/ppt_build_v1.py
# 用 LibreOffice headless 将 PPT 转为本目录同名 PDF
python code/ppt_package_v1.py
python tools/check_repository.py
```

服务器运行目录：`/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A/20260928T150000Z_ppt_redraw_v1`。上述绘图脚本固定输出本版本；新的数据或统计版本应更改 RUN/输出路径。汇总字段沿用各冻结版本，通过 `figure_manifest.tsv`、`figure_data/` 和源表清单追溯；本版没有新的检验族。
