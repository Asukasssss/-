from pathlib import Path
import json,subprocess
ROOT=Path(__file__).resolve().parents[2]
rev=subprocess.check_output(['git','rev-parse','2971100'],cwd=ROOT,text=True).strip()
A='results/PDAC/07_INTEGRATION/20260922T135400Z_report_adapter_v5/'
S='results/PDAC/06_EXTERNAL/20260922T131800Z_author_identity_v4/'
c=json.loads((ROOT/'runtime/handoff_45c6ae9/configs/BRCA_report_v1.yaml').read_text(encoding='utf8'))
c.update(cancer='PDAC',cohort='CAMP / GSE62452',input_commit=rev,protocol_commit='45c6ae9ab489a9bd4ead4e47d28185785a448156',expected=dict(metabolites=307,workpool=51,relations=357,genes=687,current_genes=250),assets=[],subtypes=[],subtype_examples=[],legacy_source_display_columns=[])
c['files']={k:A+v for k,v in dict(statistics='statistics_long.tsv',genes='candidate_genes.tsv',reader='candidate_genes_reader.tsv',relations='direct_candidate_relations.tsv',all_relations='candidate_relations.tsv',metabolites='metabolite_results.tsv',mapping='metabolite_mapping_status.tsv',rna='rna_results.tsv',rna_history='rna_history_supplement.tsv',association='association_results.tsv',conditional_association='conditional_association_results.tsv',stability='sc_source_status.tsv',cross_study='sc_cross_study.tsv',coverage='sc_gene_coverage.tsv',progress='progress_summary.tsv',gaps='module_gaps.tsv',identity='author_identity_counts.tsv',sample_audit='sample_audit_summary.tsv').items()}
cohorts=['GSE263733','GSE278688','GSE242230']
for s in cohorts:c['files'][s]=S+'sc_celltype_profiles_'+s+'.tsv'
c['columns']['relation_id']='relation_id'
c['families']=dict(metabolite_primary='primary',metabolite_sensitivity='both_preimputation_available',association_primary='primary',association_sensitivity='availability')
c['sc'].update(profile_sources=cohorts,stability_sources=['stability'],partition='primary_tumor',cohorts=cohorts,stability_top_field='top_celltype',heatmap_page_rows=26,celltypes=['Malignant (author)','Normal epithelial (author)','Ductal (unresolved)','Acinar','Fibroblast/CAF','Endothelial','Myeloid','T/NK','B/Plasma','Mast','Erythrocyte','Platelet'],labels={'Malignant (author)':'作者恶性','Normal epithelial (author)':'作者正常上皮','Ductal (unresolved)':'导管（身份未定）','Acinar':'腺泡','Fibroblast/CAF':'成纤维/CAF','Endothelial':'内皮','Myeloid':'髓系','T/NK':'T/NK','B/Plasma':'B/浆细胞','Mast':'肥大细胞','Erythrocyte':'红细胞','Platelet':'血小板'})
c['modules']=dict(rna=True,single_cell=True,subtypes=False)
c['design']=[dict(label='肿瘤 / 正常标本',value='27 / 12'),dict(label='明确代谢物与RNA配对',value=11),dict(label='关联作者病例键',value=21),dict(label='单细胞候选并集',value=687)]
c['design_source']='results/PDAC/03_PATIENT/20260922T053500Z_internal_paired51_v2/README_CN.md'
c['examples']=[dict(gene=g,relation_id=r,role=role,reason=why,limitation=lim) for g,r,role,why,lim in [
 ('ABAT','3fc531a39ec1bfa2','内部关联支持','GABA关系在主分析与可用值均过FDR；来源不一致','同队列选择；三队列来源不同，不是独立验证'),
 ('AASS','3b65c5b363143034','未支持示例','用户关注；直接生化关系但患者统计较弱','RNA P=0.0527；相关未支持；α-KG原有值仅4对'),
 ('MGLL','d47762fd25c589b8','版本边界示例','2-palmitoylglycerol关系保留；新版直接族q=0.059','不能引用旧冻结42入口的q=0.02715冒充本轮'),
 ('SLC6A19','720d1bf4e114a655','来源身份示例','色氨酸直接运输关系；GSE242230首位作者恶性','主关联q=0.0531；另两队列导管身份未定'),
 ('TDO2','385ef55bc5987f6d','稳定来源示例','三队列成纤维/CAF首位且保持率达门槛','主关联q=0.0856；来源稳定不证明作用细胞或机制')]]
c['narratives']={
 'design':['本次哪些样本能用于哪种分析？','39个跨组学标本：27肿瘤、12正常；11个作者明确配对用于代谢物与RNA，21个不同作者病例键用于肿瘤内关联。','另6个肿瘤缺病例键，未默认独立。临床身份和分装层级未再认证。'],
 'metabolites':['全量代谢物中有哪些探索线索？','307项均可评估；51项P<0.05，0项q<0.05。入口是探索性工作池。原CAMP冻结303/42另行保留，不与307/51混算。','配对精确符号秩和原10000次整对区间复用；处理尺度不是浓度倍数。'],
 'paired_counts':['配对变化方向有多少一致？','展示固定候选对应代谢物的全部11对方向计数。α-KG为6升、1降、4相等；比例分母包含相等对。','相等值可能受到作者处理影响；人数不是逐患者显著次数。完整排序与所有代谢物见附录。'],
 'sensitivity':['限制到作者原有值后，效应如何变化？','193/307项达到至少8对门槛；27项P<0.05，0项q<0.05。其余114项保留不可评估状态。','作者表可用性不是质谱检出认证；缺项不能当阴性。均值区间不等于符号秩检验反演区间。'],
 'mapping':['51项怎样形成候选范围？','357条直接关系对应250基因；358条条件关系另列。条件及历史并集687基因全部进入来源，7项无可靠关系仍保留。','图中非直接基因437个包括条件与历史补充，不全是历史独有。命名和身份限制不被映射消除。'],
 'associations':['肿瘤之间代谢物和RNA是否共同变化？','357条直接关系中354条可评估：65条P<0.05、1条q<0.05。条件关系另族350/358可评估，4条q<0.05。','21个作者病例键，未调整协变量；同队列探索。直接、条件及可用值各有独立BH族。'],
 'forest':['固定示例的相关及敏感性怎样？','GABA—ABAT主rho=0.849、q=0.0354，可用值q=0.0342。AASS—α-KG主rho=-0.223，可用值rho=-0.014。','示例含支持、未支持、身份及版本边界；不是按最小q形成的靶点排名。'],
 'rna':['AASS等候选在配对RNA中如何变化？','当前250基因中247可评估，91项P<0.05、40项q<0.05；补充437基因独立一族，429可评估、0项q<0.05。','RNA来自已处理连续微阵列，配对t检验；丰度不等于酶活性或通量。'],
 'rna_counts':['RNA的患者内方向是否一致？','全部方向来自同一11对。AASS为8对降低、3对升高，P=0.0527，q=0.1369。完整250基因结果含NA。','方向人数不能取代配对检验；作者尺度差值不称log2倍数。RNA未显著仍进入单细胞。'],
 'dot':['候选在哪些细胞类别表达？','全部687基因已按来源标签等权汇总。AASS偏CAF/内皮，TDO2为三队列CAF首位；各研究分别显示色标。','GSE242230拆分10783作者恶性与350正常上皮；另两研究导管身份未定。最高表达不证明唯一作用细胞。'],
 'evidence':['不同证据如何并列？','代谢物、内部关联、可用值及RNA各自保留原P/q。直接主分析只有GABA—ABAT过FDR，不能将其他层显著移作该关系支持。','所有列不合成评分；单细胞没有疾病差异检验P/q。'],
 'discussion':['当前可以提出哪些问题？','保留ABAT的内部联系、AASS的未支持结果、MGLL/SLC6A19的版本边界和TDO2的来源线索。完整715关系及687基因保留。','机制、临床分型和外部同关系验证不由本报告补证；停在来源。'],
 'gaps':['哪些模块目前不能展示？','两队列恶性身份、临床分型连接、作者UMAP坐标和跨研究患者独立性仍有缺项。点图与热图覆盖完整基因池。','不以细胞Classical/Basal状态代替患者分型；不凭肿瘤组织来源认定恶性，不临时新建嵌入。'],
 'stability':['来源一致和可靠定位覆盖怎样？','三研究可解释来源分别621、593、610/687；120基因三队列首位相同，72个同时保持率≥80%，当前直接池30个。','缺项显示暂不可定位；共同类别另表。不同恶性/正常/身份未定类别不强行合并。']}
c['mapping_labels']={'DIRECT_MAPPED':'直接映射','CONDITIONAL_ONLY':'仅条件性','UNRESOLVED':'暂未映射','NAMED_UNRESOLVED':'身份待核','UNKNOWN_FEATURE':'未知特征'}
for n in range(1,16):
 for ext in ['png','svg']:
  c['assets'].append(dict(role='appendix',path=S+f'figures/dotplot_{n:02d}.{ext}',source_table=S+'figures/figure_manifest.tsv',caption=f'v4全基因点图第{n}页，原字节复用。共用色标仅保留历史展示，不用于跨研究绝对表达比较；新版主报告各研究单独色标。'))
(ROOT/'configs').mkdir(exist_ok=True)
(ROOT/'configs/PDAC_report_v5.yaml').write_text(json.dumps(c,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print(rev)
