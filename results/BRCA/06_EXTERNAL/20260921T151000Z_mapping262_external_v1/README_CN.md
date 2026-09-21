# 新版262关系外部关联

## 本轮问题
为新版262条关系提供FUSCC、Tang外部患者证据。全部条目保留，不只检验CAMP显著项。

## 输入与范围
沿用明确作者样本连接：FUSCC 258例TNBC、Tang 20例当前RNA版本可连接病例；同源数据不重复算独立验证。源矩阵仍在server165。

## 实际结果
|家族|计划|可评估|复用|新算|新q<0.05|
|---|---:|---:|---:|---:|---:|
|FUSCC_EXCLUDE_WARNING262|262|126|125|1|12|
|FUSCC_PRIMARY262|262|126|125|1|13|
|Tang262|262|235|153|82|1|

FUSCC新增RNA接口拒绝连接，32个基因的请求受阻；涉及主分析64条已匹配代谢物但缺RNA的关系，标ACCESS_BLOCKED。旧缓存继续使用。Tang新版唯一q<0.05关系为MDH1—苹果酸；旧效应/P未改，PNP原线索保留。

## 新手解释
扩展名单后，旧关系的P值和效应保持不变，但新版q会随检验范围变化。q是否过0.05的改变不是生物学作用突然出现或消失。

## 限制与反证
FUSCC是阶段性可用覆盖：接口恢复后补齐受阻项，再另建完整q版本。此次BH按262个计划项、缺项内部p=1处理，公开缺项p/q仍NA；Tang按235个可评估检验校正。两队列q不能当作统一重要性分数。
作者处理代谢组包含既有填补，原始仪器缺失标志未额外核实。名字或标识匹配不等于重新确认化合物身份。LPCAT4待稳定基因身份核对。相关不能证明代谢介导、酶活或因果；Tang小队列区间不精确。Oslo真实编号连接仍未解决，不猜连接。

## 当前决定
保留未支持、缺测、受阻结果；已有ASNS/GLS/SLC6A8等线索继续保留。新增候选不能因FUSCC尚未下载RNA被淘汰。

## 下一步
与CAMP、新39功能材料回接；接口恢复时只补缺口。新版q不能覆盖旧174结果。

## 复现
独占server165运行目录，source/direct_relations.tsv放置新版映射，.running内容mapping262_external_v1。
`python3 brca_mapping262_external_v1.py <new_run_directory>`
脚本code/brca_mapping262_external_v1.py。公开标记和独立BH复核由code/brca_publish_mapping262_v1.py external生成。
