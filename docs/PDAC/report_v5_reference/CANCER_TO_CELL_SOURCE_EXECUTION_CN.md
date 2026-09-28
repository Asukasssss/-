# 各癌种执行合同 v1：CAMP发现至单细胞表达来源

这是执行交接，不是再次审核全包的要求。用户要求各癌种达到与BRCA相同的完成范围、记录结构和展示质量。先复用已完成且来源明确的计算，只补缺失步骤。禁止为追求相同阳性数量复制BRCA样本数、基因名单、结果或结论。

## 0. 开工、权限与完成范围

- 每个窗口只处理用户已分配的癌种；当前用户在该窗口的明确分配优先于可能过时的assignments.tsv，记录冲突后更新自己分配，不接管其他负责人。无分配时只询问癌种，不自动选癌种。
- 仓库 https://github.com/Asukasssss/- 。先fetch，读取AGENTS.md、docs/STAGE_STANDARD_CN.md、自己癌种的阶段索引和最近README。使用自己的clone/工作树与`analysis/<cancer>-<task>`分支。
- 本合同与展示程序目前在独立交接分支，未自动合并main；先确认能读取其提交。展示程序参照提交8952c54b61709c875ec83881b66ee43bd1de59f3。不要把BRCA整套结果复制成自己的分析。
- server165目录使用`results/collaborative/<CANCER>/<ACCOUNT>/<RUN_ID>/{private,public}`并加独占运行锁。原始矩阵、样本连接、供者明细、逐细胞坐标/矩阵、凭据只留服务器；本地和GitHub只交公开汇总、图件、代码与说明。
- 终点：样本设计→代谢物发现→直接基因映射→CAMP肿瘤内部关联＋RNA背景→全候选单细胞来源→整合→真实报告。单细胞到来源停止。
- 功能文献、DepMap、外部患者复现、CPTAC、生存、拟时序、通讯、机制干预不作为本轮前置条件；已有结果链接保留，禁止为了补齐阶段新增这些任务。
- 不统一制造配对。无可证实配对或正常RNA可走登记的替代路径；没有适合单细胞资料则阶段PARTIAL/ACCESS_BLOCKED，不宣称已完成定位。两套单细胞是优先目标，不能取得第二套时先交一套，明确复现缺口。

## 1. 目录与执行顺序

保留既有七阶段目录，不另发明冲突编号。下面是执行顺序，阶段编号是登记位置，不是先后关卡。

|执行顺序|工作|登记位置|标准化汇总文件名|
|---|---|---|---|
|0|样本与设计确认|03_PATIENT；原CAMP输入留01_CAMP|design_summary.tsv、sample_audit_summary.tsv|
|1|代谢物主分析及缺失敏感性|01_CAMP保留冻结结果；新增配对/敏感性留04_ROBUSTNESS|metabolite_results.tsv、metabolite_workpool.tsv|
|2|直接关系及历史基因并集|02_MAPPING|metabolite_mapping_status.tsv、direct_relations.tsv、conditional_relations.tsv、gene_membership.tsv|
|3|肿瘤内部关联与RNA背景|03_PATIENT；新增敏感性索引到04_ROBUSTNESS|association_results.tsv、rna_results.tsv|
|4|单细胞表达来源|06_EXTERNAL，analysis_type注明sc_expression_source|sc_feature_coverage.tsv、sc_celltype_profiles.tsv、sc_source_status.tsv、sc_cross_study.tsv|
|5|统一候选与展示|07_INTEGRATION|statistics_long.tsv、candidate_relations.tsv、candidate_genes.tsv、candidate_genes_reader.tsv|

目录均为`results/<CANCER>/<stage_id>/<RUN_ID>/`。已有历史文件无需改名或挪动；用source_manifest及适配列映射到上述逻辑表。同一统计只存一份，其他阶段用索引引用。05_FUNCTION无新增则标NOT_RUN或引用历史，不制造结果。

每阶段必须有README_CN.md、analysis_spec.json、source_manifest.tsv、validation.json、公开明细；更新`coordination/stages/<CANCER>.tsv`，沿用templates/stage_index.tsv表头，stage_id/run_id排序。完成一批就commit/push并ls-remote核对，不等整个癌种结束。

## 2. 每一步的计算与记录

### 第一步：确认实际统计单位

分别确认：肿瘤—正常同患者配对、RNA—代谢组同标本/同患者连接、重复区域/时间点、患者与标本数量、癌种病理、处理前后与队列来源。配对必须来自作者明确对照或可追溯元数据，不能按排序、数字后缀或分型猜。

私有逐样本连接留服务器；公开只记录各层样本数、匹配规则出处、冲突数/处理规则、排除原因计数。患者独立性未确认时写author_sample/source_donor_label，不能直接写独立患者。

多个队列分别分析；重复患者必须先确定一患者一单位规则或适合重复测量的模型。不能把同一患者多个区域当作独立重复。癌种没有配对时先登记设计分支，继续可做的肿瘤内部关联和单细胞。

### 第二步：全量代谢物发现

冻结原CAMP效应、P/q、分类，另列新分析版本。分析全量特征，不能先挑已显著部分才做BH。

**明确配对且沿用BRCA可比的处理值时：**令d=肿瘤−同患者正常，使用双侧Wilcoxon符号秩检验。BRCA参照实现为zero_method=wilcox、method=approx、correction=True；至少8个完整独立配对，非零差少于10标近似警示，全零差P=1。其他癌种若小样本/大量结/差值分布不适合该实现，计算前冻结精确/置换或其他合适方法并说明差异，不机械照抄。符号秩的位置解释依赖差值对称性。

同时报告平均配对差、中位差、秩二列效应、n升/n降/n相等、n升/n及n降/n。均值差的95%区间用整对病例bootstrap4000次百分位法；这是均值的点区间，不是Wilcoxon检验的反演区间。处理尺度不清时不称logFC或浓度倍数。

主分析复用作者处理值，不擅自填补。敏感性只保留作者未填补表两侧均有值的配对，仍比较相应处理值；该掩码叫“作者表可用值”，除非有证据不能叫检出限或原始检出。无掩码则NOT_EVALUABLE，不伪造零。

每个队列的主分析和敏感性分别对全部可评估P作BH，登记计划数、可评估数、排除原因。主分析P<0.05进入探索性映射清单；保留q并标P_ONLY或FDR_SUPPORTED。全量表仍完整保留。

排序固定：升高/降低分表，各按相应方向配对比例降序、P升序、稳定特征键排序；比例分母包含相等对，另列有效非零数。均值方向与多数配对方向冲突时显式标记，不隐藏。敏感性未支持不从全量表删除。

**没有明确配对：**不得生成“配对比例”。可复用合适的原CAMP非配对统计作为探索入口；需要新比较时依据实际设计预先确定Welch/秩和/含批次的模型，记录为独立版本。分组计数和效应正常报告，配对字段NA。无法排除肿瘤/正常与批次完全混杂时只保留条件性背景，不作可靠疾病效应结论。

### 第三步：直接代谢物—基因映射

按本癌种探索入口的精确特征键连接原映射；可靠旧关系复用，新增特征补人类直接反应酶、直接运输者、必要复合体成员。未知X特征保留，身份/异构体/脂质层级不明的关系放conditional，不按名字猜基因、不扩成整个通路。

每条关系记录metabolite_key、原名称、稳定化学ID、gene/稳定基因ID、role、substrate/product/transport/complex角色、反应ID、方向/可逆性、物种、证据URL/版本、实验/相似性推断、匹配精度、当前/历史/条件性状态。多文献出处另表，一条精确关系只计一次；不同测量峰即使同名也不自动合并。

gene_membership分别列当前直接池、历史保留、新增、并集。旧基因的RNA/单细胞可按稳定基因ID复用；旧关系P/q不能移给同基因的新代谢物关系。入口P显著只是工作分批，不是映射后机制成立。

### 第四步A：CAMP肿瘤内部代谢物—RNA关联

覆盖当前全部合格直接关系。只用明确连接的肿瘤，正常不混入；逐关系完整值子集记录实际n。BRCA默认至少8个独立可评估单位，常量/缺测/连接有歧义均单列状态。

参照实现为Spearman rho；双侧9999次置换，P=(至少同样极端的次数+1)/(9999+1)，4000次按独立病例整体重抽样得到rho百分位95%点区间。随机种子由固定版本和完整关系键确定，记录有效bootstrap次数；无效重抽样不填零。重复测量或交换性不成立时不能套独立样本置换，须先冻结适当单位/分层方案。

主分析与作者表可用值敏感性分别构成全量检验族，BH覆盖该队列该版本全部可评估关系，不只校正P<0.05或展示基因。缺项保留NA及原因；新增关系扩大范围时另列新版q，相同输入/方法的效应和P可复用。

保留P<0.05的探索性视图并区分q<0.05，但单细胞必须覆盖全部当前候选及历史保留基因，不能以相关过线再次截断基因池。未做协变量调整就明确unadjusted，不能继承BRCA其他模型结论。本轮不为了过线反复加协变量模型。

### 第四步B：RNA疾病背景

与相关分析分开。配对证据由RNA自己的样本表确认，不能照搬代谢组配对。作者已处理、连续且尺度适用时，BRCA参照为配对t检验：均值d、t分布95%区间、4000整对bootstrap辅助区间，n升/降/相等及比例；至少8对且差值变异可评估。RNA尺度不明不得称log2FC。

如果是原始RNA-seq整数计数，采用适合计数的edgeR/DESeq2设计（例如patient+condition），记录过滤、归一化、设计和对比；不能把原始计数直接配对t检验。非配对连续表达可用预先指定Welch/limma等，批次可辨识时纳入模型。无正常RNA标NOT_EVALUABLE，继续其他部分。

RNA独立BH族为本队列当前候选池全部可评估基因；保留P<0.05工作视图及全部P/q。历史保留基因的旧结果另列版本，不与当前q混算。RNA升降不是抑制/激活靶点建议。

### 第五步：单细胞只看基因表达来源

优先选本癌种有作者注释、表达矩阵、可追溯供者/样本元数据的两套肿瘤研究；检视患者重叠。正常参考可选，不能算第二套肿瘤验证。没有第二套先完成第一套，状态明确。

1. 覆盖全部当前基因＋历史保留基因，不仅看精选或相关显著者；用稳定基因ID核对探针/别名。一对多或身份冲突单项暂停。
2. 保留作者细胞注释和来源。宽类别映射表逐项保存原label→统一label，无法对应标Unresolved，不静默丢细胞。不把非恶性上皮统称癌细胞；其他癌种不强制都有乳腺癌类别。
3. 有原始counts时采用log1p(count/每细胞全基因库大小×10000)，库大小不能只按候选基因求和。已有合适作者标准化矩阵可复用并记录尺度，不能重复log/归一化。整合表达/残差层不自动用于丰度及检出比例。
4. 先每个供者×细胞类计算细胞平均表达、counts>0检出比例和细胞数，再在供者间等权平均，避免细胞多的患者独占结果。只有源标签时注明source_donor_label，不能冒称确认独立供者。
5. 参照覆盖下限：每个供者×类至少20细胞，每类至少3个有效供者；不足记录NOT_EVALUABLE。参数可因实际数据调整，但必须在结果前登记理由和版本，不以阳性数调参。
6. 可解释排名要求至少两个可评估类别且最高平均检出比例≥1%。保留机械top/runner作底表，但不合格时易读表显示“暂不可定位”，不显示“两研究来源一致”。覆盖条件与低表达条件分列。
7. 按供者联合重抽样1000次，保留缺类情况，计算原top保持率、有效重抽样次数；最高并列按并列数量分摊。报告top/runner平均差及同时覆盖两类供者中的方向比例。保持率不是患者一致百分比，也不是显著性P。
8. 两研究分别报告全部可评估类别top、共同类别top、覆盖数与状态。两研究原top相同且各保持率≥80%仅标描述性稳定；不同平台不直接比较绝对颜色/均值。缺项不是反向证据。
9. 有可靠临床分型可分层重复来源描述：分别保存各型自身类别、各型共同类别、等权分型汇总三种口径。至少满足对应覆盖条件。临床亚型不等于单细胞状态标签；不照搬ER/HER2/TNBC。没有分型标NOT_RUN。本轮不新增分型差异显著性检验。
10. UMAP优先复用作者或已有坐标：表达与细胞类别图使用同一细胞集合、坐标和方向；标处理尺度、细胞数、注释来源、抽样规则。无坐标可先交点图/热图，UMAP缺项明确；如另建一次探索嵌入须单独登记输入、QC、特征、参数和种子，不冒称作者坐标。不为每个基因单独重算嵌入。

不追加拟时序、细胞通讯、通路机制、推断代谢通量或“来源细胞即唯一作用细胞”的结论。

## 3. 所有窗口通用的字段合同

数值表前缀沿用templates/statistical_result.tsv，不能任意换列顺序：

`cancer, cohort, stage_id, run_id, analysis_version, analysis_type, metabolite_key, metabolite_name, gene, unit, n, n_reference, effect_type, effect, ci_lower, ci_upper, p_value, q_value, test_family, family_n_evaluable, status, reason, source_id`

追加字段而非覆盖：relation_id、稳定基因/特征ID、direction、升/降/相等人数与比例、planned_n、缺失数、seed、bootstrap有效数、reused_from、input_sha256、evidence_url、annotation_source。精确统计键至少包括队列×真实测量特征×基因×分析类型×版本；独立reaction来源不能制造重复检验。

单细胞没有检验的P/q填NA：追加study、partition、celltype、n_donors、n_cells、mean_expression、mean_detection、expression_scale、source_status、top_raw、top_display、runner、top_frequency、valid_bootstrap_n、shared_category_n。公开表保留汇总，不导出供者ID或每供者数值。

状态统一DONE/PARTIAL/NOT_RUN/NOT_EVALUABLE/ACCESS_BLOCKED/NEEDS_REVIEW；reason必填缺项原因。真实0与NA区分。统计“未显著”、技术“未测到”、工作“未运行”不可混写。

analysis_spec必须包含：输入版本/哈希、单位、纳排规则、设计、尺度/变换、缺失处理、方法/对比、检验族计划/可评估数、随机种子、软件版本、复用范围、探索性选择边界。README固定“问题→输入范围→实际结果→新手解释→限制/反证→当前决定→下一步→复现”。

## 4. 整合与推理记录

- candidate_relations按精确关系保留代谢物主/敏感性、肿瘤相关主/敏感性、基因RNA、来源状态、历史外部/功能链接；不将某基因最佳代谢物换来换去后比较模型。
- candidate_genes保留当前/历史并集；易读表显示P支持与FDR支持分层，不以P/q大小合成靶点评分。所有原表在附录完整保留。
- INTERPRETATION_CN.md分开“直接观察到的数值”“可提出的解释”“不能推断的内容”“下一步最小问题”。每项引用关系ID/图源，标专门讨论或通用保留说明，不把通用文字冒称深入研究。
- 不强制所有层都显著。表达最高不证明作用最强，RNA相关不证明反应方向/通量，某癌种显著另一个不显著不证明癌种特异。

## 5. 与BRCA一致的真实展示交付

主报告建议12–16页，框架固定，缺失模块仍占位说明，不伪造图凑页：

|顺序|页面内容|
|---|---|
|1–3|目标与边界、样本设计、各阶段覆盖/计数|
|4–6|全量代谢物效应与P/q层级、配对方向比例或非配对替代、映射/历史并集|
|7–9|关系相关效应及区间、敏感性、RNA差异与方向|
|10–12|两研究来源点图/热图、来源稳定性/缺项、已有分型或缺失说明|
|13–15|同坐标表达/类型UMAP示例、证据矩阵、研究问题和限制|

每页写“问题—图—实际结果—解释—限制”。精选仅为展示，在display_selection.tsv先登记完整关系ID及理由，包含支持/未支持/缺项示例；不按最小q临时凑名单，完整候选始终在附录。

必须交付：PDF＋Markdown、可筛选/冻结表头Excel、每张新图PNG＋PDF/SVG＋图源TSV＋中文图注、全基因来源热图/点图、全基因UMAP图册（有坐标时）、figure_manifest/source_manifest/checksums、代码/配置/依赖/验证、可解压完整ZIP。原栅格图复用不伪造矢量。

颜色/版式沿用BRCA白底、中文解释、缺项灰色；同一细胞类别跨本癌种图一致。不同表达尺度分开图例；热图行缩放须标只用于同基因跨细胞类别比较。UMAP不画未经统计支持的“显著富集”圈。

复用code/camp_results_report.py，参照docs/CAMP_RESULTS_REPORT_USAGE_CN.md。先为本癌种建立字段适配和配置，完整更换cancer/cohort/design/files/families/sc/expected/examples/narratives及输入完整SHA。当前程序已完成BRCA真数据与合成缺项测试，不保证任意癌种表一键兼容；不适配时修正读取适配，不改科学数值。生成报告前冻结自己的输入提交，不使用运行时HEAD/main。

```
python code/camp_results_report.py --config configs/<CANCER>_report_v1.yaml --validate-only
python code/camp_results_report.py --config configs/<CANCER>_report_v1.yaml --out results/<CANCER>/07_INTEGRATION/<NEW_RUN_ID>
```

## 6. 验收和给用户的最终回复

逐项填写templates/cancer_to_source_acceptance.tsv，不把缺项改成DONE。最低核对：键唯一；配对方向计数和=n；每族BH范围固定；q不跨版本挪用；全候选未漏；来源状态屏蔽机械top；分型口径分开；历史值不变；逐页无乱码/截断/缺图；ZIP能解压、哈希一致；GitHub不含受限数据。

每癌种公开一张progress_summary.tsv：步骤、计划数、实际可评估数、P<0.05数、q<0.05数、status、缺口原因、路径。无P/q的来源阶段填NA，不写0。

最终必须报告：做到哪一步、每步真实数量/单位、缺失模块、未完成原因、重点结果及推理限制、主PDF/Excel/ZIP链接、输入完整SHA、远程交付提交、是否合并、复现命令。真正停止在单细胞来源和整合，不自行继续机制挖掘。仅有流程文档或代码不等于完成，必须给本癌种实际生成结果。
