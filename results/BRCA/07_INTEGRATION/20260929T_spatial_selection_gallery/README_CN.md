# LYPLA1 四癌空间图选片库
本轮仅整理已有公开衍生图，不重算统计，不根据结果筛除图。
BRCA 82 幅、COAD/CRC 31 幅、PDAC 4 幅、PRAD 12 幅，共 129 幅。数量包含汇总、分辨率版本及排版版本，不是独立切片数。
COAD 文件夹包含 CRC 研究，不将直肠切片称为纯 COAD。PDAC GeoMx 为区域比较；Moncada H&E 尚未完成像素级配准。GSE235315 检查样本无 LYPLA1 条目，不伪造表达图。
输入路径与每幅 SHA256 见 source_manifest.json。原图、矩阵与历史统计不修改。
预览：D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/LYPLA1_spatial_selection_all/index.html
复现：python code/spatial_selection_gallery.py
支持按癌种筛选、搜索、放大、勾选、复制编号。当前未完成浏览器交互验收：浏览器工具启动失败。文件存在、图像解码与哈希另行验证。
