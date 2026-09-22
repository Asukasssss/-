# PDAC 51 项配对探索工作池映射 v2

## 本轮问题
按用户指定“配对发现→直接映射→CAMP内部关联/配对RNA→三套单细胞→外部验证”推进；不以旧42条限制新工作池。

## 输入与范围
复用307项、11对配对分析，51项原始P<0.05、均q≥0.05。依据作者原化学注释，结合KEGG/ChEBI/Rhea与reviewed人类UniProt反应注释；保留原字段与来源链接。新增18个KEGG化合物查询，全部缓存哈希核对。

## 实际结果
51项：30直接关系、12仅条件、3身份暂挂、6范围内未找到可靠精确关系。715条去重关系=357直接+358条件，共538个候选基因。mapping_evidence.tsv保留逐反应证据，planned_relations.tsv为唯一待测关系；条件支持行可能有多个出处，不按出处重复检验。

## 新手解释
这是可检验的关系台账，不是715条已验证机制。SAH/alpha-ketoglutarate涉及蛋白或核酸修饰的共底物/产物注释放入条件层；复合体成员单列，不称独立催化酶。51条只达到未校正P门槛，都是探索性。

## 限制/反证
化学身份未重新实验认证；Pipecolate的D/L、1-stearoylglycerol原药物混合物ID和游离C-glycosyltryptophan身份限制保留。未找到关系只指当前检索范围，非不存在。注释不是全面逐关系功能审查。

## 当前决定
固定357/358两个关联族及全部538基因RNA/单细胞池，再运行统计；不依据新关联结果回改规则。详情见docs/PDAC/PAIRED51_INTERNAL_V2_LOCK.md。

## 下一步
21个作者明确且互不重复配对编号的肿瘤做内部关联；同11对做RNA背景；全部基因接三队列来源。

## 复现命令
`python code/pdac/build_mapping_paired51_v2.py`（需原版本来源缓存；输出目录必须新建）。来源及参数详见mapping_sources.tsv、source_manifest.tsv、analysis_spec.json。
