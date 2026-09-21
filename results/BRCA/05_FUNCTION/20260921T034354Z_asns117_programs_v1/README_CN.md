# BRCA：117候选回接Asns干预与三项转录程序

## 本轮问题
哪些原117候选响应Asns敲低？预先固定的氨基酸饥饿应答、氨基酸跨膜运输及EMT相关表达是否发生方向一致的变化？本轮是GSE104966同研究公开RNA结果再利用，不是新的独立实验验证。

## 输入与范围
复用7f77c0b发布的四份完整edgeR结果和全部基因覆盖表，不重拟合患者或RNA计数，不改变原49,246项联合BH q。原发瘤来源与肺转移来源均为小鼠4T1-T经过体内生长、分离、6TG筛选及再培养的细胞；每个来源每组4个作者重复，两条构建共享相应对照。独立小鼠身份未新增核实，不推断配对。

人—鼠同源使用Ensembl110，要求唯一one-to-one；多个人源返回项通过稳定ID lookup的精确正式名称消歧，不按大小写猜测。该身份规则细化发生在初次映射发现同义名冲突后，未改变DE或程序检验参数；见mapping_resolution_note.json。

程序使用MSigDB2025.1.Mm原生小鼠集合，Asns从集合和全部已检验背景中删除。limma3.58.1 cameraPR对有方向sqrt(QLF F)做竞争性秩检验，固定基因间相关0.01，未从样本估计。3集合×4比较统一BH12，最少10个覆盖成员。详细成员及原q全部保留。该程序检验属于同研究、后选择探索。

## 实际结果
117个全部保留，111个唯一一对一映射；ADA2、AHCY、CS、PNP、PRODH、SLC6A8未满足该规则，不能因此判为没有功能或绝对没有同源。原发瘤来源87个可评估、24个低表达未检验、6个映射未决；肺来源88个可评估、23个低表达未检验、6个映射未决。共468行完整候选×来源×构建记录，234行两构建汇总。

|两构建响应分类|原发瘤来源|肺转移来源|
|---|---:|---:|
|均显著且同向下降|15|1|
|均显著且同向上升|8|0|
|同向、仅一条构建显著|18|3|
|同向、均未过原联合q阈值|39|57|
|点估计方向相反|7|27|
|未检验或映射未决|30|29|

原发瘤23个共同显著者包括干预对象ASNS自身，因此下游候选为22个：
- 上升8个：CEPT1、FDFT1、GOT1、MAT2A、SLC6A6、SOAT1、SPHK1、SQLE。
- 下降14个：AOX1、ARG2、CMPK1、CPT1A、CROT、GAMT、GSR、LDHB、NAGK、PANK2、PPA2、RFK、SORD、UPP1。

肺来源仅ASNS自身符合两构建共同显著同向规则；这不是其他候选“没有响应”的等效性结论。

|程序|原发瘤sh1 q|原发瘤sh2 q|肺sh1 q|肺sh2 q|方向说明|
|---|---:|---:|---:|---:|---|
|氨基酸饥饿应答|0.3107|0.2389|0.6278|0.6231|均未获本轮支持|
|氨基酸跨膜运输|0.002977|0.04947|0.3107|0.9057|原发瘤两构建偏上；肺未支持|
|EMT相关表达|0.01414|0.05604|1.572e-9|0.01907|原发瘤偏上但仅一构建过阈值；肺两构建偏下|

图：program_response.png、all117_response_heatmap.png。色彩为描述性log2FC（程序图为集合基因中位数），星号分别对应BH12与原49,246项联合q。程序图星号不是对中位数本身做的检验。灰色为不可评估，非零效应。

## 新手解释
“共同显著同向”表示两种降低Asns的方法都伴随该基因RNA向同一方向变化；共享对照，不能算两次独立研究。22个是Asns干预响应候选，不是22个基因各自获得干预验证。

运输相关基因的表达整体偏上，提示可以关注营养运输相关转录调整，但本轮没有直接测运输速度或营养摄取。EMT相关表达偏下不能单独证明转移受到抑制。饥饿应答未显著不能证明细胞没有应激。

## 限制/反证
仅RNA，未测蛋白、酶活、谷氨酰胺丰度或代谢通量；不能连接成ASNS—GLS功能补偿链。没有来源×干预交互检验，不能由两个来源显著基因数或方向差异直接宣称来源特异作用、耐药。程序固定相关0.01为假设，集合大小和重叠见明细；竞争性检验比较集合相对其余已检验基因的排序。表达不显著不等于变化很小，未新增等效检验或区间估计。

此前Gls/Gls2/Glul没有原联合q支持补偿性上调的结论不变。本轮统计复制逐值核对，最大误差低于1e-12；BH12另行复算一致。没有重新运行原edgeR模型或原49,246项BH。完整历史比较列逐字符串保留，见integration_validation.json。

## 当前决定
保存全117响应注释及有限程序分析，不按本轮q重新筛掉原患者/功能候选。6个同源未决单独保留。22个响应候选用于提出下游转录问题，不启动22条独立大课题。GSE104967零对照与GLS直接干预资源问题仍未解决。

## 下一步
优先取得设计明确的GLS自身干预资料；若无就保持未完成，不用营养剥夺代替。现有运输程序可用于与适用干预资料比较，暂不反复增加本研究通路或调整阈值。EMT属于原研究已有主题，不宣称首次发现。人类患者负关联和小鼠干预响应继续作为不同层证据。

## 复现命令
服务器运行目录：/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A/20260921T034354Z_asns117_programs_v1

冻结候选清单genes117.json与冻结程序成员frozen_program_sets.json已随汇总交付；用于复现时放入运行目录source/。保留服务器source缓存及上一轮public输入，在该目录运行：
```sh
python3 brca_asns117_map_v1.py --root "$PWD"
Rscript brca_asns_programs_v1.R "$PWD"
python3 brca_asns117_validate_plot_v1.py --root "$PWD"
```
R命令采用服务器已安装limma环境；详见代码。取回public汇总后在仓库运行：
```sh
python code/brca_finalize_asns117_v1.py
python tools/check_repository.py
```

主要文件：all117_comparison_asns_response_appended.tsv为历史比较表追加7列；all117_four_contrast_response.tsv为四比较原统计复用表（logFC=RNA log2FC，p_value=原P，original_q_global4=原联合q，status/reason=可评估性）；不是新增117项检验。program_results.tsv采用统一统计表前缀，新q仅属于BH12。source_manifest.tsv及code_manifest.tsv记录输入与脚本哈希，源矩阵和API缓存留server165。

来源：Ensembl https://e110.rest.ensembl.org ；MSigDB https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Mm/ ；GSE104966 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104966 。
