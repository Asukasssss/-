"""One-time BRCA configuration authoring; not part of the generic renderer."""
from pathlib import Path
import json,subprocess
R=Path(__file__).resolve().parents[1]
P='results/BRCA/07_INTEGRATION/20260922T112000Z_sop_record_v1/'
M='results/BRCA/04_ROBUSTNESS/20260921T112611Z_metabolite318_paired_v1/'
S='results/BRCA/06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/'
N='results/BRCA/06_EXTERNAL/20260922T023000Z_sc39_v1/'
L='results/BRCA/06_EXTERNAL/20260919T140105Z_scRNA117_v1/'
A='results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/'
files={'statistics':P+'statistics_long.tsv','genes':P+'candidate_genes156_history_preserved.tsv','reader':P+'candidate_genes156_reader.tsv','relations':P+'candidate_relations262_history_preserved.tsv','metabolites':M+'paired_metabolite318.tsv','mapping':'results/BRCA/02_MAPPING/20260921T124136Z_mapping190_v2/metabolite190_mapping_status.tsv','rna':'results/BRCA/03_PATIENT/20260921T150000Z_mapping262_patient_v1/paired_RNA150.tsv','association':'results/BRCA/03_PATIENT/20260921T150000Z_mapping262_patient_v1/CAMP_relations262.tsv','subtypes':S+'genes156_subtype_comparison.tsv','subtype_profiles':S+'subtype_celltype_profiles.tsv','history_inventory':P+'historical_artifact_inventory.tsv'}
examples=[
 ('LYPLA1','KEGG:C04102|LYPLA1','患者联系','主关联与可用值分析支持；恶性上皮背景','精确脂质外部关系未评估，不是游离GPC'),
 ('ASNS','KEGG:C00064|ASNS','患者联系','RNA升高及谷氨酰胺负关联','ER+浆细胞排名与类别覆盖限制；非通量'),
 ('GPI','KEGG:C00668|GPI','反证示例','内部正关联但配对RNA未支持','Tang已评估未支持正关联'),
 ('GPCPD1','KEGG:C00670|GPCPD1','反证示例','内部负关联及髓系背景','外部支持不足；RNA配对未支持'),
 ('NNMT','KEGG:C02918|NNMT','背景敏感性','正关联与成纤维背景','组成调整减弱；RNA配对降低'),
 ('LYPLA2','KEGG:C04102|LYPLA2','RNA/来源示例','RNA升高及恶性上皮背景','此具体关系不因RNA显著而升级'),
 ('ABHD12','KEGG:C00219|ABHD12','RNA/来源示例','RNA升高及髓系来源','不以名义相关替代功能验证'),
 ('LPCAT1','KEGG:C04102|LPCAT1','来源差异示例','RNA升高','两研究来源不同；分型第一第二名接近'),
 ('SLC7A11','KEGG:C00491|SLC7A11','RNA/来源示例','多数配对RNA升高','两研究最高类别不同；关系不强行解释'),
 ('PYCR1','KEGG:C00148|PYCR1','RNA/来源示例','RNA升高及上皮背景','已有CAF功能文献不能被最高类别覆盖'),
 ('AASS','KEGG:C00047|AASS','患者联系','RNA降低、正关联、成纤维背景','敏感性q跨阈值但rho不变；外部不一致'),
 ('ETNK1','KEGG:C00189|ETNK1','患者联系','主分析及可用值支持','来源不同；作用依赖背景，非既定干预方向'),
]
sc=[]
for name in ['Wu2021','Pal2021_reprocessed']:
 for folder in [L,N]:
  alias=name+('_old' if folder==L else '_new')
  files[alias]=folder+name+'_celltype_profiles.tsv';sc.append(alias)
files['stability_old']=A+'all117_source_stability.tsv';files['stability_new']=N+'new39_source_stability.tsv'
config=dict(cancer='BRCA',cohort='CAMP / Terunuma',input_commit=subprocess.check_output(['git','rev-parse','8663004'],cwd=R,text=True).strip(),protocol_commit=subprocess.check_output(['git','rev-parse','6dd0a73'],cwd=R,text=True).strip(),expected=dict(metabolites=318,workpool=190,relations=262,genes=156,current_genes=150),files=files,
 columns=dict(gene='gene',relation_id='relation_key',metabolite_id='metabolite_key',metabolite_name='metabolite_name',effect='effect',p='p_value',q='q_value',status='status',n='n',analysis='analysis_type'),
 families=dict(metabolite_primary='paired_processed',metabolite_sensitivity='paired_author_available',association_primary='processed_Spearman262',association_sensitivity='available_Spearman262'),
 sc=dict(profile_sources=sc,stability_sources=['stability_old','stability_new'],partition='ALL',cohorts=['Wu2021','Pal2021_reprocessed'],celltypes=['B_cells','Endothelial','Fibroblasts','Malignant_epithelial','Mast_cells','Myeloid','Nonmalignant_epithelial','Perivascular','Plasma_cells','T_cells'],labels={'B_cells':'B细胞','Endothelial':'内皮','Fibroblasts':'成纤维','Malignant_epithelial':'恶性上皮','Mast_cells':'肥大细胞','Myeloid':'髓系','Nonmalignant_epithelial':'非恶性上皮','Perivascular':'血管周','Plasma_cells':'浆细胞','T_cells':'T细胞'},partition_field='partition',cohort_field='cohort',type_field='celltype',detection_field='mean_detection_fraction',stability_top_field='top_lineage',stability_max_detection='max_detection',heatmap_page_rows=26),
 subtypes=['ER+','HER2+','TNBC'],subtype_examples=['ASNS','LPCAT1','SHMT2','LYPLA1','ABHD12','ENPP2'],
 examples=[dict(gene=g,relation_id=k,role=role,reason=reason,limitation=limit) for g,k,role,reason,limit in examples],
 design=[dict(label='当前跨组学标本',value=108),dict(label='明确肿瘤—正常配对',value=45),dict(label='肿瘤内部关联作者病例',value=60),dict(label='组织标签冲突标本',value=1)],
 design_source='results/BRCA/03_PATIENT/20260921T102429Z_camp_sample_identity_v1/README_CN.md',
 mapping_labels={'DIRECT_MAPPED':'直接映射','CONDITIONAL_ONLY':'仅条件性','NAMED_UNRESOLVED':'已命名未解决','UNKNOWN_FEATURE':'未知特征'},
 font_candidates=['C:/Windows/Fonts/msyh.ttc','C:/Windows/Fonts/simhei.ttf','/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'],
 modules=dict(rna=True,subtypes=True,single_cell=True),
 assets=[dict(role='main_umap',path='results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/LYPLA1_expression_celltype_pair.png',caption='LYPLA1表达与同坐标细胞类型参照。使用已有栅格图与原色标，不重新计算坐标或恢复矢量。',source_table='results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/gene_page_index.tsv'),
 *[dict(role='appendix',path='results/BRCA/06_EXTERNAL/'+run+'/'+file,caption=label,source_table='results/BRCA/06_EXTERNAL/'+run+'/'+table) for run,file,table,label in [
 ('20260922T082327Z_sc156_umap_v1','all156_expression_UMAP.pdf','gene_scales_and_pages.tsv','全156基因表达图；原全细胞坐标'),
 ('20260922T085000Z_sc156_umap_pairs_v1','all156_expression_celltype_pairs.pdf','gene_page_index.tsv','全156基因表达与类型对照'),
 ('20260922T092000Z_epithelial_umap_v1','all156_epithelial_triptychs.pdf','gene_page_index.tsv','全156基因上皮三联图；此前单独计算的探索嵌入，非全细胞坐标子集')]]],
 narratives={
 'design':['哪些测量可以对应？','作者对应表支持配对设计；配对数和肿瘤内部关联病例数是不同统计单位。','相同作者标本编号不证明同一分装；组织冲突标本已在新关联中排除，冲突本身仍待解释。'],
 'metabolites':['肿瘤相对自身正常组织，哪些代谢物变化？','全量展示效应与P/q层级；探索入口依据P，完整表保留所有项目。','均值差沿用作者处理尺度，不是浓度倍数；同队列探索，不是外部验证。'],
 'paired_counts':['示例代谢物有多少配对方向一致？','升高、降低和相等分别计数；分母取每项实际可用配对数。','每对只表示方向，不代表对每位患者分别做显著性检验。'],
 'sensitivity':['可用值子集是否改变效应？','主分析和敏感性并排；样本量与缺项保留，不能选择更有利的版本。','作者data可用值不是新增核验的仪器检测掩码；图为汇总效应，不是患者散点。'],
 'mapping':['代谢物如何形成候选池？','直接映射、条件性、未解决与未知特征全部有去处；当前与历史基因分开。','直接关系包含酶、转运体及必要复合体成员，不等于功能已验证。'],
 'associations':['同一肿瘤中的代谢物与RNA是否共同变化？','展示全部计划关系，未测和待核对项目保留。','rho是关联强度，不是倍数；候选来自同队列，不能称独立验证。'],
 'forest':['不同证据类型的示例关系表现怎样？','展示已保存的rho与点95%区间，并列主分析、敏感性。','12个示例不是前十名或最终靶点；RNA/来源示例的关系不因被展示而升级。'],
 'rna':['同患者肿瘤中的RNA是否发生变化？','展示全部当前候选及P/q分层；P入口仅用于展示。','使用已保存配对t结果与作者表达尺度，不把RNA差异当酶活。'],
 'rna_counts':['RNA方向是否在患者间一致？','选择固定示例，展示升高、降低和相等配对数量。','RNA方向一致不代表靶点验证；RNA不显著也不阻止单细胞定位。'],
 'dot':['这些基因在哪些细胞背景表达？','颜色表示该研究供者等权平均log1p表达，点面积表示平均检出比例；不可评估标叉。','各研究单独色标，不跨平台比较绝对颜色；表达最高不是唯一作用细胞。'],
 'subtypes':['分型的来源排名是否一致？','自身可用类别与共同类别分别展示；SHMT2为ER+髓系、另两型恶性上皮。','等权汇总最高不等于各亚型都最高；ASNS类别覆盖不齐；LPCAT1第一第二名接近。'],
 'umap':['高表达区域与哪些细胞群重合？','现成表达图与相同坐标的类型图一起阅读；完整候选图册见附录。','UMAP不是富集检验；图像不能替代供者级来源表。'],
 'evidence':['证据怎样并列而不机械取交集？','代谢物、关联、敏感性、RNA分别保留P/q；来源和外部限制另列。','颜色仅标统计层级，不构成统一分数；缺测不能记成0。'],
 'discussion':['这份结果最适合讨论什么？','LYPLA1、ASNS患者线索较清楚；AASS、ETNK1也有具体联系；RNA/来源示例同样保留。','GPI的Tang未支持、GPCPD1外部不足、NNMT背景敏感及研究覆盖限制均不能隐藏。']})
config['sc']['labels'].update(T_NK_cells='T/NK细胞',Unresolved_stromal='未定基质')
config['legacy_source_display_columns']=['Wu_top','Pal_top']
config['files'].update(normal_reference_profiles=L+'Reed2024_celltype_profiles.tsv',normal_reference_coverage=L+'Reed2024_feature_coverage.tsv',old_source_comparison=L+'gene117_cell_source_comparison.tsv',new_source_comparison=N+'new39_cross_study_comparison.tsv',subtype_coverage=S+'subtype_celltype_coverage.tsv',subtype_within_celltype=S+'within_celltype_subtype_comparison.tsv')
conclusions={
 'design':'108个标本中，45对支持配对比较；组织冲突排除后的肿瘤内关联为60个作者病例。',
 'metabolites':'318项全部可评估：190项P<0.05，其中176项q<0.05。190项用于探索映射，不能全部称FDR支持。',
 'paired_counts':'谷氨酰胺39/45对降低；1-palmitoyl-GPC 36/45对升高；两者构成方向不同的示例。',
 'sensitivity':'可用值子集288项可评估、121项q<0.05；两种分析均显著且均值同向的为93项。',
 'mapping':'190项中93项有直接关系、13项仅条件性、31项已命名未解决、53项未知；得到262关系和150当前基因。',
 'associations':'262条计划关系中255条可评估：69条P<0.05、14条q<0.05；7条缺测或身份待核。',
 'forest':'LYPLA1—1-palmitoyl-GPC主rho=0.482、可用值rho=0.518；AASS的rho不变而q跨阈值，不能称效应消失。',
 'rna':'150个当前基因中148个可评估：85个P<0.05、79个q<0.05。RNA差异不作为单细胞入场门槛。',
 'rna_counts':'SLC7A11为44/45对升高，LYPLA1为34/45对升高，AASS为37/45对降低；方向人数不等于功能验证。',
 'dot':'LYPLA1/LYPLA2偏恶性上皮，ABHD12偏髓系，NNMT偏成纤维；LPCAT1两研究最高类别不同，不能强指定作用细胞。',
}
for key,value in conclusions.items():config['narratives'][key][1]=value
(R/'configs').mkdir(exist_ok=True)
(R/'configs/BRCA_report_v1.yaml').write_bytes((json.dumps(config,ensure_ascii=False,indent=2)+'\n').encode())
print('BRCA configuration written')
