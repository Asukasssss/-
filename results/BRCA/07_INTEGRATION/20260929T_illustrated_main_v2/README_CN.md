# CAMP 泛癌发现与 LYPLA1：可编辑组会汇报

## 本轮问题
将三个代表案例之前的成果收拢成正式 PPTX，统一版式并提高可编辑性。

## 输入与范围
采用后续已确认版本：开场五页重建为原生文字、连接线与独立科研插图；队列概况、同向代谢物和基因（含患者一致性）、代谢物–LYPLA1桥接、组织配对、四癌全细胞/单核、选定BRCA/PRAD空间、近期文献、简洁证据汇总。插图来源与哈希见illustration_manifest.json；科学结果来源沿用前版。

## 实际结果
共18页，见slide_outline.tsv。没有纳入ASNS、SLC6A6、UCKL1三个案例。原始科学结果和阈值保持不变。原生对象计数及WPS打开/导出结果见validation.json。

## 新手解释
插图是独立可裁剪、可替换的位图，内部像素不作为矢量形状编辑。标题、正文、数字、图例、箭头、点阵和患者连线可在WPS/PowerPoint选择、改字或改样式；UMAP点云、病理照片、论文截图为独立可替换图片，没有整页截图拼接。矢量图形是原生可编辑形状，不是连接Excel数据的原生统计图表；改数据应从原分析代码重新绘制。

## 限制/反证
跨癌总览含六癌，LYPLA1专题限定BRCA/COAD/PDAC/PRAD。四癌组织平均上调，PDAC未显著；不能称四癌均显著。患者比例不等于每位患者显著。空间定位不等于癌细胞特异性或因果证据。两个PRAD切片来自同一患者。文献仅复用已核对页面，不代表本次再次系统检索。部分密集统计页保留小字；可在PPT中按汇报需求拆页。WPS已实际打开导出；未在Microsoft PowerPoint单独渲染。

## 当前决定
交付可编辑PPTX、WPS导出PDF、18页预览与来源说明。统一中文字体显式设置东亚字体，避免WPS自动替换为宋体。

## 下一步
组会审阅内容与讲述节奏；确认后将三个代表案例放入后续章节。

## 复现命令
python code/ppt_illustrated_main_v2.py
powershell -File code/ppt_illustrated_main_export_v2.ps1
python code/ppt_illustrated_main_package_v2.py
