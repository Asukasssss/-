# CAMP 完整组会汇报：25页

## 本轮问题
将用户指定的代表案例v3接入插图精修版。
## 输入与范围
原18页插图精修版，加CAMP_representative_cases_v3中按网页顺序的7页。19–20为ASNS，21–23为SLC6A6，24–25为UCKL1。无案例空间页。
## 实际结果
25页PPTX及WPS实际导出PDF。前18页XML逐页完全一致；新增页面统一页眉、中文字体与页码。来源及哈希见source_manifest.json。
## 新手解释
文字、反应图与患者连线是可编辑对象。UMAP点云与图像为独立可替换图片，不是整页截图。
## 限制/反证
不重算统计、不改变证据措辞。生化注释不等于本队列机制证明；相关性不代表因果。
## 当前决定
采用用户明确指定v3七页，接在LYPLA1证据汇总之后。
## 下一步
按组会时间审阅整套讲述节奏。
## 复现命令
python code/ppt_complete_cases_v3.py
使用WPS打开生成PPTX，导出25页PNG到rendered及PDF。
python code/ppt_complete_cases_package_v3.py
