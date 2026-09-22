# CAMP癌种内发现至单细胞表达来源：可复用执行规范 v1

版本：2026-09-22。范围：从CAMP输入到候选整合，单细胞仅做表达来源描述。本文是新癌种启动规范与交接模板，不是宣称所有癌种已计算，也不是可直接运行所有历史脚本的总入口。

## 0．主线、边界与版本

**样本身份核查 → 全量代谢物配对分析及缺失敏感性 → P<0.05探索工作池 → 直接生化映射 → 肿瘤内部代谢物—RNA关联＋配对RNA背景 → 全候选单细胞表达来源 → 关系级与基因级整合。**

外部患者验证与功能深入不作为本主线的开工前置；已有外部和功能证据原样接回。基础功能文献可并行，但不扩展为本规范的必做步骤。单细胞不追加细胞通讯、拟时序、RNA velocity、调控网络、通路机制或正式亚型差异检验。上皮专用UMAP属于可选展示，不是完成标准。

适用癌种按仓库固定顺序：BRCA、COAD、GBM、PDAC、PRAD、ccRCC。相同的是问题、记录结构、证据边界；配对设计、平台尺度、实际n及必要模型按队列确定。

历史数据冻结：原CAMP效应/P/q不覆盖；相同输入、样本、键和方法的结果优先复用；新候选加入时旧效应/P可复用，但新的检验集合另算并另存q。没有新计算时不得把旧分析称为本轮新增验证。

### 0.1 六步与仓库七阶段的对应

|本流程|仓库stage_id|说明|
|---|---|---|
|身份核查、全量发现|01_CAMP；必要敏感性04_ROBUSTNESS|已有BRCA配对补算位于04，保留原路径|
|直接生化映射|02_MAPPING|映射和统计分开|
|内部相关、RNA差异|03_PATIENT|敏感性可索引04，数值不重复存成两份证据|
|单细胞来源|06_EXTERNAL|属于外部细胞表达背景，不等于代谢关系外部复现|
|最终比较表|07_INTEGRATION|来源和数值版本逐列可追溯|

不重编号原STAGE_STANDARD_CN.md，不修改其他癌种结果。本文及模板作为共同规范补充，通过PR交接；各癌种执行时冻结自己的analysis_spec，而不是覆盖公共模板。

### 0.2 三个层级不能混淆

- 代谢特征：队列/平台中的实际峰或特征键。名称相同但峰不同，不能自动合并。
- 直接关系：`cohort × metabolite_key × stable_gene_id × mapping_version`。多篇来源支持同一关系只增加证据记录，不增加统计检验数。
- 基因：去重的稳定基因ID。RNA差异和单细胞按基因；相关分析按关系。

全量代谢物表、当前关系池、当前基因池、历史基因并集分别保存。不能用一张仅含显著项的表替代全部结果。

## 1．运行前固定输入、参数与记录

### 1.1 最小输入

1. 代谢物作者处理矩阵、特征注释与单位/变换说明。
2. 原有值或缺失掩码及其准确含义；没有时记录不可用。
3. RNA矩阵与基因/探针注释、处理尺度说明。
4. 患者、标本、组织、RNA与代谢组的作者连接表及证据来源。
5. 相关原CAMP结果，用于保留历史效应/P/q，不作为覆盖新全量检验范围的依据。

每个文件登记 `source_id,path_or_url,release,download_date,sha256,table_or_sheet,unit,transform,missing_encoding,imputation,access_scope`。先识别NA、0、LOD、负值分别代表什么；不得把真正零值改成缺失，或把填充值当原测得值。

### 1.2 推荐的新队列默认参数

以下是项目工作默认值，不是普遍的生物统计充分性标准。必须在看结果前写入analysis_spec；调整必须留理由和新版本。

|事项|默认值/规则|
|---|---|
|双侧显著性展示|P<0.05；BH调整P<0.05另标FDR_SUPPORTED|
|探索入口|配对代谢物主分析P<0.05；不先筛原非配对显著项|
|最低可计算配对/相关n|8；低于8保留描述与NOT_EVALUABLE；n=8也不代表充分功效|
|均值差/相关CI|4000次整患者/整对重采样，百分位95%区间；逐项区间，不是同时区间|
|Spearman P|9999次置换＋1校正；全部可评估关系使用相同预定规则|
|单细胞覆盖|每供者×类别≥20细胞，每研究×类别≥3来源供者标签|
|单细胞排名最低信号|至少两类可评估，至少一类平均检测比例≥1%；否则不强排第一|
|单细胞排名稳定性|1000次按供者标签重采样，另存有效次数；80%只是描述性标志|
|随机数|固定master_seed；SHA256(版本+队列+分析族+稳定键)派生逐项种子|

写入软件及版本、实际P算法、零差处理、ties处理、连续性校正、bootstrap种子、CI对象、检验族planned/evaluable数。不要依赖会随软件版本变化的`method='auto'`。

### 1.3 每批固定交付

`README_CN.md`、完整TSV、`analysis_spec.json`、`source_manifest.tsv`、`validation.json`、`checksums.tsv`、代码与环境版本、阶段索引。README按“问题→输入→实际结果→新手解释→限制→当前决定→下一步→复现”排序。

公开统计表以前缀开始：

`cancer,cohort,stage_id,run_id,analysis_version,analysis_type,metabolite_key,metabolite_name,gene,unit,n,n_reference,effect_type,effect,ci_lower,ci_upper,p_value,q_value,test_family,family_n_evaluable,status,reason,source_id`

新增字段接在后面。状态使用DONE、PARTIAL、NOT_RUN、NOT_EVALUABLE、ACCESS_BLOCKED、NEEDS_REVIEW。缺测为NA并给reason；0只能用于真实数值或真实计数。

## 2．步骤一：核实样本身份、真实配对和组学连接

### 2.1 逐项确认

- 患者、标本、切块/分装、实验测量ID分别建列，不能把它们都叫sample。
- 肿瘤与正常由作者明确病例表/可靠元数据连接；不靠行顺序、数字后缀、分型或表达相似性猜配对。
- RNA与代谢组确认连接层级：同一分装、同一标本不同分装、同患者不同组织块，或未知。
- 正常组织类别明确：癌旁、远端正常、健康供者等不能无说明合并。
- 去除真正技术重复或按预定规则合并；多个肿瘤区域不能当多个独立患者。优先按预定质量规则选单一代表，或用患者聚类/层级设计；不能挑最显著区域。
- 组织标签冲突单项暂挂，说明影响的配对和相关分析；不把整癌种阻塞。
- 同一患者RNA有配对，不代表其代谢物必然有配对；分别生成两种配对表。

### 2.2 输出如何记录

服务器private保存：`sample_identity_audit_private.tsv`、`metabolite_pairs_private.tsv`、`RNA_pairs_private.tsv`、`tumor_multiomics_map_private.tsv`。

必要字段：`patient_id,specimen_id,tissue_label_author,tissue_label_verified,metabolomics_id,RNA_id,pair_id,link_level,link_source,duplicate_group,decision,reason`。配对表每行一对，分别列tumor/normal测量ID。相关连接表每行一个独立分析单位。

GitHub只上传`sample_audit_summary.tsv`：`audit_item,n_checked,n_pass,n_excluded,n_unresolved,status,reason,source_id`，以及脱离个人/标本ID的决策说明。模板含私有表空表头，可上传；填入真实ID的表不可上传。

### 2.3 必须通过的局部检查

唯一键无重复；一对中同一患者且组织类型不同；每个ID在矩阵存在；one-to-one连接预期一致；排除数与总数相加一致；两组学链接不因merge重复而扩行。

**完成标准：**有可追溯连接的分析单位可启动；不能确认的单位单独记录。没有可靠配对时不执行“配对”模型，按预先指定的非配对/混合设计建独立版本；单细胞不必等待不存在的正常RNA。

## 3．步骤二：全量代谢物配对发现

### 3.1 固定范围与尺度

对该癌种输入中全部保留特征一次性分析，不仅分析原CAMP显著项。继承作者处理尺度并明确单位，不无理由再log、再填补、再批次校正。若原尺度无法支持比较，先记录并固定必要预处理的新版本。

对第i位患者和代谢物m：`d_i = tumor_i - normal_i`。只用两侧均有有效值的整对；不可分别删两列NA后按剩余行号配对。

### 3.2 描述量与方向比例

每个特征保存：总合格配对数、该特征实际配对数、升高数、降低数、相等数、均值差、中位差。检查 `n_up+n_down+n_equal=n_pairs_used`。

- `up_fraction=n_up/n_pairs_used`，`down_fraction=n_down/n_pairs_used`，默认包含持平对作为分母。
- 可另存非零差方向比例，但列名注明分母为`n_nonzero`，不可与上项混用。
- 排序视图可按`max(up_fraction,down_fraction)`降序，再按实际n降序、P升序；原全表不删除。
- 平均方向、median方向、多数对方向可能不同，各自记录并标`direction_discordant`，不能强行改成一个方向。

### 3.3 配对检验与效应

主要检验：双侧Wilcoxon signed-rank。先计算d，按数据已知精度处理数值误差；若无明确精度说明，不凭结果挑舍入位数。零差采用`wilcox`（检验时去掉，描述计数仍保留），ties用平均秩。

效应单独保存：

`mean_delta=mean(d)`；`median_delta=median(d)`；`rank_biserial=(W_positive-W_negative)/(W_positive+W_negative)`。

rank-biserial按非零绝对差的平均秩计算，不是Hedges g。主CI对应均值差，用4000次整对bootstrap；若对median或rank-biserial另给CI，应新列明确对象，不能共用同一区间。

**新队列P算法推荐冻结为：**全部零差时项目约定P=1、rank-biserial=0并标记；否则n_nonzero≤16枚举符号翻转，17–29做99999次符号翻转Monte Carlo，≥30用含ties修正与连续性修正的正态近似。符号翻转用`abs(sum(sign(d)*rank(abs(d))))`为双侧统计量；全枚举取极端比例，Monte Carlo用(b+1)/(B+1)。枚举/置换仍需差值符号可交换等假设，并不自动解决严重非对称性。

**BRCA历史不同点：**现存代谢物配对版本全部用正态近似`wilcoxon(...zero_method='wilcox',correction=True,alternative='two-sided',method='approx')`，最低8对，小非零n另有警示。本文新队列推荐的P算法不追溯替换BRCA历史P/q。若严格复现BRCA，应冻结其SciPy版本与原参数；新版SciPy可能把同类参数名称改成`asymptotic`，不能只复制字符串。

Wilcoxon对位置变化的解释依赖差值分布的对称性等条件；不要把它写成任何分布下均检验“均值等于0”。若设计不适用，预先指定其他模型并另留版本，而不是结果出来后轮换检验直到显著。[方法定义：SciPy Wilcoxon](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html)

### 3.4 缺失/填补敏感性

主分析复用作者处理矩阵。敏感性限定该特征在肿瘤、正常两侧均有作者原有值的配对，使用相同处理尺度上的数值重新计算。它主要检查填补样本的影响，不是新的独立队列。

若作者`data`表只表示“有限值可用”，未经确认不能称为原始检测掩码；字段写`author_available_value_mask`。无真实掩码时标NOT_EVALUABLE，不声称完成了缺失敏感性。

分别记录两族P/q、n、效应和区间，再用完整metabolite_key连接，记录方向一致、效应差和覆盖下降。不要仅按两边q都过线定义稳健。

### 3.5 输出与完成标准

- `metabolite_paired.tsv`：主分析全特征，包括不可评估。
- `metabolite_available_sensitivity.tsv`：敏感性全特征，包括不可评估。
- `metabolite_workpool.tsv`：主分析P<0.05的工作视图；保留P、q、效应、方向人数、缺失标签。
- `metabolite_direction_ranking.tsv`：仅展示排序，不增加检验。

附加字段至少：`n_pairs_total,n_pairs_used,n_nonzero,n_up,n_down,n_equal,up_fraction,down_fraction,mean_delta,median_delta,rank_biserial,n_both_available,missing_rate_tumor,missing_rate_normal,mask_semantics,p_method,zero_method,tie_method,continuity_correction,B,seed,bootstrap_valid,ci_method,ci_estimand`。

完成时核对全特征覆盖、升降持平计数、方向与差值一致、BH范围、至少一个独立公式/工具复算。P<0.05表示探索入口，不是全部通过FDR。

## 4．统一的多重校正与版本规则

每个队列分别固定检验族；不把不同问题或不同癌种的P混为一个未经规划的集合。

|test_family|计划集合|q如何计算|
|---|---|---|
|METAB_PAIRED_PRIMARY|全代谢特征|全可评估P一次BH|
|METAB_PAIRED_AVAILABLE|同一全特征的可用值敏感性|独立一次BH|
|REL_TUMOR_PRIMARY|当前全部去重直接关系|全可评估P一次BH|
|REL_TUMOR_AVAILABLE|同一关系池敏感性|独立一次BH|
|RNA_PAIRED_CURRENT|当前全部去重基因|全可评估P一次BH|
|SC_SOURCE_DESCRIPTIVE|无假设检验|P/q=NA，不套用RNA q|

默认`family_n_planned`保留计划总数，`family_n_evaluable`为合法P数；不可测项P/q=NA。若历史版本使用缺项按P=1占位，应原样保存并注明，不因结果改变规则。

BH：排序后的P为p_(i)，`q_(i)=min(1,min_{j>=i}(m*p_(j)/j))`。这里字段q指BH调整P，不是Storey q-value或“此基因为假的概率”。使用`multipletests(...,method='fdr_bh')`或`p.adjust(...,method='BH')`并显式记录方法。标准BH保证依赖相应独立/正依赖条件，不能声称任意依赖下均成立。[R官方校正说明](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)

标签：FDR_SUPPORTED（q<.05）；NOMINAL_EXPLORATORY（P<.05且q≥.05）；NOT_SUPPORTED_THIS_TEST（已测但未满足上述标准）；NOT_EVALUABLE（未测/不可算）。后两类不是同一概念。

只对显著子集再BH是错误的。加入新关系后可以保留旧P/效应，重新校正新全族，分别存`q_original,q_current,test_family,analysis_version`。不同检验族q不能当统一重要性分数。

由于候选入口来自同一CAMP的差异分析，后续同队列相关和多次探索都要标`same_cohort_postselection`；各节点BH不等于控制整个项目所有选择的总体错误率。CI为逐项区间，不能把事后挑出的最强关系CI称为选择校正区间。

## 5．步骤三：从工作代谢物映射直接人类基因

### 5.1 身份先于关系

保留作者特征ID、原名称、标准名称、ChEBI/HMDB/其他稳定标识、构型、脂质链长/不饱和度/位置层级、身份可信度与来源。只允许有记录的名称规范化，不能模糊匹配后悄悄当作精确匹配。

例如游离GPC与1-palmitoyl-GPC不是同一个化合物；D/L或1-/2-位置未明确时不猜。X未知峰保留差异，不能据名称猜基因。

### 5.2 合格关系类型

1. 直接反应酶：该明确代谢物是反应底物或产物。
2. 直接转运者：有针对该底物的人类运输注释/证据。
3. 必要复合体成员：与直接反应的必要复合体有明确关系，角色注明complex_member，不叫独立催化酶。

数据库/原始论文至少记录来源ID、版本或访问日、物种、反应式或运输对象、证据级别。优先可靠人类条目和明确反应；Rhea/UniProt等用于追溯，不把数据库列出统一写成人类直接实验证明。`By similarity`另列INFERRED。

通路邻居、转录调控、间接信号关系不混入DIRECT；底物层级或身份未解决为CONDITIONAL/UNRESOLVED。基因别名先回到稳定ID，不能仅做大小写转换。

### 5.3 表结构与去重

- `metabolite_identity.tsv`：一行一个原特征。
- `mapping_evidence.tsv`：一行一条关系的一个来源，允许同关系多出处。
- `direct_relations.tsv`：一行一个合格去重关系，`relation_id`唯一。
- `conditional_relations.tsv`：条件关系及缺口。
- `gene_pool_current.tsv`与`gene_pool_history_union.tsv`：分别为当前池和历史并集。

关系字段至少：`relation_id,metabolite_key,stable_gene_id,gene,relation_type,metabolite_role,reaction_id,reaction_text,organism,evidence_id,source_url,evidence_level,mapping_status,identity_status,compound_specificity,added_or_reused,mapping_version`。

同名基因的不同代谢物关系不合并统计。旧基因的单细胞/功能背景可回接，新关系不能继承旧关系的相关P/q。身份变化、基因更名必须保留旧ID与新ID及证据。

完成标准：工作池每个代谢物均有状态；DIRECT去重；当前基因集合等于DIRECT基因去重；历史并集无丢失；条件关系不误入主统计族。映射完整不等于生化关系穷尽。

## 6．步骤四A：CAMP肿瘤内部代谢物—RNA关联

### 6.1 分析单位

仅取组织标签可靠且两组学正确连接的肿瘤。不是把肿瘤和正常混在一起相关，也不是默认计算`ΔRNA`与`Δ代谢物`。

每关系按实际有效两列取交集n；有正常对照不是纳入条件，因此关联肿瘤数可多于配对数。多标本来自同患者时按已冻结的代表样本/聚类策略处理，不把区域当独立患者。

### 6.2 计算

Spearman：两列分别转平均秩，`rho=Pearson(rank(M),rank(RNA))`。n≥8，且两列各至少两个不同值才可算；常数列记NOT_EVALUABLE，不能填rho=0。

默认P用患者单位置换：固定M，打乱RNA在独立患者之间的对应；保留列内分布和ties。9999次，双侧`P=(1+count(abs(rho_perm)>=abs(rho_obs)))/(B+1)`。分层或重复设计不能无条件全局打乱，应采用预先制定的受限置换/模型。置换只检验相应可交换性下的无关联，不提供因果结论。

bootstrap按患者整行重采样，同次抽样一起带走M和RNA，4000次，取rho的2.5/97.5分位。重采样产生常数列时该次无效；记录valid/total，默认低于90%有效则CI记NA并给原因，rho/P可保留。

近似P与置换P不混为同一字段。若为复用历史统计采用近似P，明确p_method。SciPy提醒小样本Spearman近似P精度有限，适合考虑置换。[SciPy Spearman](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html)

9999次Monte Carlo最小非零P为0.0001，q在阈值附近有模拟误差；记录B及seed，不换种子追求显著。如计划采用更多置换，应在查看结果前冻结或对整个指定族统一新版本，保留原值。

### 6.3 缺失敏感性与解释

主分析作者处理值；敏感性只保留代谢物原有值且RNA有效的同一批肿瘤，重算完整关系池，各自BH。相同关系比较n、rho、CI、P/q；不因敏感性不显著自动删除。

正相关不等于合成，负相关不等于消耗；RNA不等于蛋白、活性或通量。rho接近0且CI很宽应叫不精确，不能简单叫无联系。

### 6.4 输出

`tumor_association.tsv`、`tumor_association_available.tsv`、`relation_coverage.tsv`。附加列：`relation_id,stable_gene_id,n_tumor_linked,n_complete,n_unique_patients,unit_verification,n_metabolite_available,rho,ci_method,bootstrap_valid,p_method,B,seed,statistics_reused,reused_from,mask_semantics,same_cohort_postselection`。

CI主前缀effect_type=Spearman_rho，effect=rho；P/q分别保留。未确认患者独立时n写实际标本数，并以unit和unit_verification明确，不擅自写患者数。

完成标准：每条DIRECT关系均有一行和状态、精确身份连接、全族BH、缺项原因。可计算关系先交付，未配套RNA的一项不阻塞全部。

## 7．步骤四B：候选基因配对RNA背景

### 7.1 与相关分析独立

该问题是同患者肿瘤相对正常的RNA变化，不要求代谢物在这些配对都可用。当前全部去重基因一次分析；未显著仍进入单细胞。

### 7.2 按RNA输入选择模型

**作者已规范处理、尺度明确的连续表达（BRCA当前情形）：**令d=RNA_T−RNA_N。默认配对t检验，`t=mean(d)/(sd(d)/sqrt(n))`，自由度n−1；均值差95%CI用t分布。附加4000次整对bootstrap均值CI作为预定敏感性，不能冒称两者检验对象完全不同。

记录均值差、中位差、升降持平数与比例。n<8或差值零方差时主t检验不算，记原因，不返回无限t的“完美显著”；描述仍保留。差值分布严重异常时预先设定替代方法或敏感性，不逐基因看P再挑模型。

**原始RNA-seq计数：**不把原counts直接输入配对t。使用完整计数矩阵的低表达过滤、TMM规范化及负二项GLM/QL等适用流程，设计至少包含患者block与组织，如`~patient+tissue`；确认设计满秩。报告log2FC及相应模型统计，过滤、离散度估计和归一化基于适用完整背景，不只拿候选几行估计。候选BH范围与全转录组q分别命名，不能混用。此分支需另存平台适配脚本和模型配置。[edgeR官方用户指南](https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)

**没有正常RNA/没有配对证据：**该项NOT_EVALUABLE；如适合做非配对RNA比较，另建明确设计版本，不伪称配对；肿瘤内部相关与单细胞继续。

基因多探针映射：稳定探针ID和注释先确认，按预定且不依赖效应/P的规则折叠；不挑P最小探针。存在旧别名歧义的单项NEEDS_REVIEW。

### 7.3 输出

`paired_RNA.tsv`（当前全基因）、`RNA_identity_coverage.tsv`、`RNA_P005_view.tsv`（展示视图，不替代全表）。

附加列：`stable_gene_id,expression_platform,expression_scale,n_pairs_total,n_pairs_used,n_up,n_down,n_equal,up_fraction,mean_delta,median_delta,model_formula,p_method,ci_method,bootstrap_mean_lower,bootstrap_mean_upper,probe_mapping_status,current_pool,history_only`。

P<.05视图保留q分层；RNA没有差异不取消对应代谢关系或单细胞分析。历史池基因继续保留旧结果，是否纳入本次RNA检验族必须预先声明；不悄悄把“当前150”改为“历史156”重算q。

## 8．步骤五：单细胞仅到表达来源

### 8.1 数据选择与身份

优先两项独立研究，核实患者来源不重复；有可读表达矩阵、供者/来源标签、作者细胞注释；记录肿瘤/正常、原发/转移、治疗、临床亚型、平台和数据处理版本。一套就绪先算，不等第二套；第二套缺失不是第一套阴性。

来源标签不自动等于独立患者；重复切块按患者合并或记录限制。作者注释、后续研究重注释、本轮推断标签分开。跨研究先固定粗类别映射表，不看结果再改分类使其一致。Cancer Basal SC之类细胞状态不等于患者临床TNBC。

候选全量进入，不能只取患者关联显著基因。稳定基因ID唯一匹配；同名多个ID不能取表达最高者。缺测基因保留NOT_EVALUABLE；未测不是零表达。

### 8.2 规范化与两层汇总

有原始UMI计数时，可继承当前项目方法：`z_cg=log(1+10000*count_cg/library_c)`，library用保留原始基因全集之和，不能只用候选基因和。原始零值检测比例以`count>0`定义。

若仅有作者已规范化矩阵，读取作者变换说明并原样利用合适尺度；不再重复log。若不能确定零值对应检测与否，检测比例记不可评估，不在缩放/中心化矩阵上用>0当作检出。不同规范化尺度的研究独立解释，不能把绝对数值硬拼到同一统计检验。

第一层（服务器private）：每个`donor × celltype × gene`计算细胞数、mean(z)、检测比例。存在重复分装时先按患者和既定组织范围合并，不能让同患者多切块变多个独立供者；不同治疗/时点不随意混合。

第二层（公开汇总）：每供者该类≥20细胞才合格；至少3合格来源标签才展示该类别均值。对合格供者均值等权平均：`mean_d(mean_cells(z))`；检测比例也先供者内算，再供者等权。保存median、q25、q75和实际n。

因此不能把全体细胞直接平均后称为“供者等权”。细胞图可以展示所有细胞，但统计汇总权重仍按供者。

### 8.3 表达位置与稳定性

每研究每基因记录最高、次高类别及差距；不足两类、最大平均检出<1%或无信号，标无法解释排名，仍保留表达表。相等/数值近似相等记录tie，不靠字母顺序声称明确第一。

可选但建议的描述性稳定性：从原来源供者标签集合整供者有放回抽样1000次；同供者的所有细胞类型一起带走，保持其内部关联。每次重算等权均值和最高类别；记录有效抽样数与原最高类别保持频率。重采样重复标签是bootstrap权重，不是新增独立供者。若该次可评估类别不足两类则无效。

最高类别并列时分数平均分给并列类别，固定浮点相等容差。可记录逐供者留一排序、共同供者中top−runner差值与方向比例；这些是选定最高类别后的描述，不称独立显著检验。

两研究比较同时给：

- `same_top_all_categories`：各研究自身全部可评估类别的最高类别是否相同。
- `same_top_shared_categories`：只在两研究共同可评估类别中重排后是否相同。
- `same_top_both_bootstrap_ge080`：前一种相同且两边保持率≥80%，只是描述性标签。

保持率80%不等于80%患者均符合，更不等于通过功能验证。两研究注释/类别范围不同导致的排名变化不能直接称生物学差异。

### 8.4 分型和展示的停止边界

作者分型连接可靠时可给“亚型×类别”的上述描述，覆盖阈值不降低来凑图。低覆盖为灰/NA；不按表达反推亚型。此处不增加亚型显著性P/q，不宣称某亚型驱动总体结果。

必备图：供者等权表达点图（颜色均值、大小检出比例）、全候选细胞类别热图。可选图：现成UMAP表达＋同坐标类型参照；按分型拆图使用相同基因色标和坐标。

若另算上皮等子集UMAP，明确新坐标、输入范围、选特征、PCA、随机种子和软件参数；不能冒充作者原图或为了漂亮岛屿反复调参。原标签是否保留与是否重新聚类分开记录。UMAP上细胞数、密度和距离不作统计富集或因果证据。

**至此停止单细胞分析。**不自动加CellChat、拟时序、SCENIC或机制通路；不把“最高表达细胞”当唯一作用细胞，也不能确定对应代谢物就在该细胞生成或消耗。

### 8.5 输出

私有：`sc_donor_profiles_private.tsv`、原矩阵、逐细胞标签/坐标。

公开：`sc_dataset_registry.tsv`、`sc_annotation_map.tsv`、`sc_gene_coverage.tsv`、`sc_celltype_profiles.tsv`、`sc_source_stability.tsv`、`sc_cross_study.tsv`、图和逐图说明。

公开统计附加字段：`stable_gene_id,partition,celltype,n_source_labels_total,n_source_labels_eligible,n_cells_total,n_cells_eligible,mean_detection_fraction,median_donor_mean,q25,q75,annotation_origin,normalization,top_celltype,runner_celltype,top_gap,bootstrap_top_frequency,bootstrap_valid,paired_source_labels,paired_positive_fraction`。

P/q全部NA，test_family=SC_SOURCE_DESCRIPTIVE；CI如有须写来源为供者重采样、属于描述而非选择校正。

## 9．步骤六：整合成可讨论的候选，不再重复筛选

### 9.1 两张主表＋长表

1. `candidate_relations_integrated.tsv`：一行一条关系。原代谢物配对P/q、相关P/q、RNA P/q分别加前缀；患者数分别保留，不能统一填一个n。
2. `candidate_genes_integrated.tsv`：一行一个历史并集基因。当前/历史状态、对应关系数量/完整关系ID、RNA背景、单细胞背景、已有功能/外部信息指向。
3. 队列/细胞类型/分型结果保留规范长表，主表通过ID和路径连接，不把多条关系压成一个未经说明的best。

同一基因在不同模型的best可能是不同代谢物；比较敏感性必须锁定完整关系ID。RNA和sc可按稳定基因ID回接，相关效应只按完整关系键回接。合并必须一对一/多对一验证，不允许默默扩行。

### 9.2 标签使用多个独立维度

`metabolite_evidence`、`patient_association_evidence`、`RNA_background`、`cell_background`、`data_limitations`、`next_action`分开，不编一个不透明总分。

一个对象可同时是“患者联系有支持”和“髓系背景”；微环境不是比癌细胞低一级。内部RNA弱但有合适功能依据的对象继续保留。缺资料与已测未支持分开。

建议工作去向：优先讨论；保留线索；需补身份/数据；历史参照。未预设具体效应阈值时不把R-strong等主观标签写成算法筛选结论。

### 9.3 完成条件

输入代谢物均有去处；当前全部关系均有患者分析状态；当前/历史全部基因均保留；单细胞每研究/基因有覆盖状态；原列无修改；所有图能回到实际表和代码；限制与反证不被简表空白隐藏。

完整比较表足以进入导师讨论，不要求所有证据层显著，不要求全部未知峰完成身份鉴定，不要求机制闭环。

## 10．目录、验收与GitHub交付

服务器：`results/collaborative/<CANCER>/<ACCOUNT>/<UTC_RUN_ID>/{source,private,public}`，独占`.running`；完成后本运行锁改为完成标记，不动他人目录。

仓库：`code/`保存脚本，`results/<CANCER>/<stage_id>/<RUN_ID>/`保存public产物，`coordination/stages/<CANCER>.tsv`登记。跨步骤共用数值只存一次，索引引用即可。

每阶段：完成数值→必要局部核查→公开/私有分离→SHA256→git commit→push→读取远程ref确认→报告链接。上传分支不等于合并main。患者ID、完整样本连接、逐患者测量、逐细胞矩阵/坐标、凭据不进公开GitHub或本地交接包。

验收清单：

- 输入版本、单位、规范化与掩码语义明确；键无重复；join不扩行。
- 计划项=可评估+各类不可评估；n与方向计数一致。
- 原统计原样保留；新q族明确；bootstrap/置换种子明确。
- 单细胞归一化分母为完整基因库；供者等权；灰色表示缺测不表示0。
- 历史池并集和列值一致；按效应而非P选图的偏倚边界明确。
- 脚本中不能残留BRCA专用45、60、117、156、262的断言用于其他癌种；这些数应由该癌种冻结输入决定。

## 11．BRCA真实实现与复用边界

|环节|可查看的仓库脚本|不能直接照搬的部分|
|---|---|---|
|身份核查|code/brca_camp_sample_identity_v1.py|作者表结构、具体冲突、45/60等计数|
|配对代谢物|code/brca_metabolite318_paired_v1.py|318项、45对、旧近似P算法、作者sheet|
|关系与RNA|code/brca_mapping262_patient_v1.py|262/150范围、BRCA RNA尺度、阻塞别名、60/45样本|
|相关方法参考|code/reference/camp_per_cancer_repro.py|调用时实际参数/版本；不是所有队列通用已验证入口|
|单细胞汇总|code/brca_sc117_profile_v1.py、code/brca_sc39_v1.py|117/39硬编码、Wu/Pal类别映射与路径|
|来源稳定性|code/brca_all117_stability_v2.py|117硬编码、既有标签与覆盖规则|
|分型展示|code/brca_sc156_subtypes_v1.py|ER+/HER2+/TNBC是BRCA标签；不是跨癌固定类别|

BRCA已完成的P/q和结论不因本文新队列默认值改变。另一癌种可复用算法与表结构，必须做薄适配层并记录实际设计；没有通用执行器就不能声称“只改癌种名一键完成”。

## 12．给另一个执行账号的启动指令

> 阅读本SOP、AGENTS.md、STAGE_STANDARD_CN.md与癌种分工。只处理被分配癌种。先列输入文件、实际样本设计和缺失语义，冻结analysis_spec；使用templates/camp_discovery_source_v1的空表和参数模板。完成身份连接后按六步运行，能算的一批立即交付，缺项保留状态。候选入口P<0.05，q独立保留；RNA和相关不作为进入单细胞的硬门槛。单细胞仅做到表达来源。已有结果按精确键复用，不修改历史P/q。不运行BRCA脚本里的硬编码路径/计数，不启动其他癌种或额外机制任务。每阶段提交公开结果与脚本到自己的分支并核对远程。

## 方法资料

- [SciPy Wilcoxon官方文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html)：零差、ties、近似与置换方法。
- [SciPy Spearman官方文档](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.spearmanr.html)：秩相关及小样本P的说明。
- [R p.adjust官方文档](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)：BH/BY等调整。
- [edgeR官方用户指南](https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)：计数模型与配对设计。

这些文档用于方法依据；实际复现必须以冻结软件版本和参数为准，不能直接依赖持续更新网页的默认值。
