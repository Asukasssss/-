"""Generate empty interchange schemas, never patient data or analysis results."""
from pathlib import Path
import csv,json,hashlib
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'templates/camp_discovery_source_v1';OUT.mkdir(parents=True,exist_ok=True)
PREFIX=(ROOT/'templates/statistical_result.tsv').read_text(encoding='utf-8').strip().split('\t')
tables={}
def add(name,grain,scope,cols,stat=False):tables[name]=dict(grain=grain,scope=scope,columns=(PREFIX if stat else [])+cols.split())
add('sample_identity_audit_private','one author specimen linkage','SERVER_PRIVATE','cancer cohort patient_id specimen_id tissue_label_author tissue_label_verified metabolomics_id RNA_id pair_id link_level link_source duplicate_group decision reason')
add('metabolite_pairs_private','one explicit patient pair','SERVER_PRIVATE','cancer cohort pair_id patient_id tumor_metabolomics_id normal_metabolomics_id pair_source status reason')
add('RNA_pairs_private','one explicit RNA patient pair','SERVER_PRIVATE','cancer cohort pair_id patient_id tumor_RNA_id normal_RNA_id pair_source status reason')
add('tumor_multiomics_map_private','one independent tumor analysis unit','SERVER_PRIVATE','cancer cohort analysis_unit_id patient_id specimen_id metabolomics_id RNA_id link_level unit_verification source_id status reason')
add('sample_audit_summary','one audit item, no identifiers','PUBLIC_AGGREGATE','cancer cohort audit_item n_checked n_pass n_excluded n_unresolved status reason source_id')
add('source_manifest','one input file','PUBLIC_PROVENANCE','source_id path_or_url release download_date sha256 table_or_sheet unit transform missing_encoding imputation access_scope')
meta='n_pairs_total n_pairs_used n_nonzero n_up n_down n_equal up_fraction down_fraction mean_delta median_delta rank_biserial n_both_available missing_rate_tumor missing_rate_normal mask_semantics p_method zero_method tie_method continuity_correction B seed bootstrap_valid ci_method ci_estimand family_n_planned direction_discordant'
for name in ['metabolite_paired','metabolite_available_sensitivity','metabolite_workpool','metabolite_direction_ranking']:add(name,'one metabolite feature within a named analysis family','PUBLIC_AGGREGATE',meta+' evidence_tier workpool_entry',True)
add('metabolite_identity','one source feature','PUBLIC_ANNOTATION','cancer cohort metabolite_key author_feature_id author_name canonical_name stable_compound_id stereochemistry lipid_resolution identity_status source_id identity_version reason')
rel='cancer cohort relation_id metabolite_key stable_gene_id gene relation_type metabolite_role reaction_id reaction_text organism evidence_id source_url evidence_level mapping_status identity_status compound_specificity added_or_reused mapping_version reason'
add('mapping_evidence','one evidence source per relation','PUBLIC_ANNOTATION',rel)
add('direct_relations','one unique direct relation','PUBLIC_ANNOTATION',rel)
add('conditional_relations','one unresolved conditional relation','PUBLIC_ANNOTATION',rel)
for name in ['gene_pool_current','gene_pool_history_union']:add(name,'one stable gene identity','PUBLIC_ANNOTATION','cancer stable_gene_id gene in_current_pool in_history_pool history_only mapping_version relation_ids source_id status reason')
assoc='relation_id stable_gene_id n_tumor_linked n_complete n_unique_patients unit_verification n_metabolite_available rho ci_method bootstrap_valid p_method B seed statistics_reused reused_from mask_semantics same_cohort_postselection family_n_planned evidence_tier'
for name in ['tumor_association','tumor_association_available']:add(name,'one exact metabolite-gene relation in one cohort/family','PUBLIC_AGGREGATE',assoc,True)
add('relation_coverage','one relation and cohort','PUBLIC_AGGREGATE','cancer cohort relation_id metabolite_key stable_gene_id gene metabolite_identity_status RNA_identity_status n_linked n_complete status reason source_id')
rna='stable_gene_id expression_platform expression_scale n_pairs_total n_pairs_used n_up n_down n_equal up_fraction mean_delta median_delta model_formula p_method ci_method bootstrap_mean_lower bootstrap_mean_upper probe_mapping_status current_pool history_only family_n_planned evidence_tier'
for name in ['paired_RNA','RNA_P005_view']:add(name,'one gene within RNA contrast/family','PUBLIC_AGGREGATE',rna,True)
add('RNA_identity_coverage','one requested gene and platform','PUBLIC_ANNOTATION','cancer cohort stable_gene_id gene measured_feature_id mapping_policy mapping_version status reason source_id')
add('sc_dataset_registry','one scRNA study release','PUBLIC_PROVENANCE','cancer cohort study_accession release matrix_layer expression_scale annotation_origin donor_identity_level disease_scope treatment_available subtype_available independence_group source_id status reason')
add('sc_annotation_map','one author class to coarse class','PUBLIC_ANNOTATION','cancer cohort author_celltype coarse_celltype annotation_origin mapping_evidence version status reason')
add('sc_gene_coverage','one requested gene and study','PUBLIC_ANNOTATION','cancer cohort stable_gene_id gene source_feature_id source_symbol match_method status reason source_id')
add('sc_donor_profiles_private','one donor-celltype-gene-partition','SERVER_PRIVATE','cancer cohort donor_id partition celltype stable_gene_id gene n_cells mean_log1p10k detection_fraction raw_count all_gene_library_sum source_id')
add('sc_celltype_profiles','one gene-celltype-study-partition','PUBLIC_AGGREGATE','stable_gene_id partition celltype n_source_labels_total n_source_labels_eligible n_cells_total n_cells_eligible mean_detection_fraction median_donor_mean q25 q75 annotation_origin normalization',True)
add('sc_source_stability','one gene-study-partition','PUBLIC_AGGREGATE','cancer cohort gene stable_gene_id partition n_evaluable_celltypes top_celltype runner_celltype top_gap max_detection tie_status bootstrap_top_frequency bootstrap_valid bootstrap_total loo_top_frequency paired_source_labels paired_positive_fraction paired_mean_delta seed status reason source_id')
add('sc_cross_study','one gene and explicit study pair','PUBLIC_AGGREGATE','cancer gene stable_gene_id study_A study_B top_A top_B n_common_categories same_top_all_categories top_A_shared top_B_shared same_top_shared_categories bootstrap_top_A bootstrap_top_B same_top_both_bootstrap_ge080 status reason')
add('candidate_relations_integrated','one relation with explicit source versions','PUBLIC_AGGREGATE','cancer cohort relation_id metabolite_key stable_gene_id gene mapping_version metabolite_effect metabolite_p metabolite_q metabolite_n metabolite_source_run association_rho association_p association_q association_n association_source_run RNA_effect RNA_p RNA_q RNA_n RNA_source_run sc_source_run metabolite_evidence patient_association_evidence RNA_background cell_background data_limitations next_action status reason')
add('candidate_genes_integrated','one gene in history union','PUBLIC_AGGREGATE','cancer stable_gene_id gene in_current_pool history_only relation_ids n_relations RNA_background RNA_source_run cell_background sc_source_run functional_source_ref external_source_ref data_limitations next_action status reason')
catalog=[]
for name,t in tables.items():
 assert len(t['columns'])==len(set(t['columns'])),name
 with (OUT/(name+'.tsv')).open('w',encoding='utf-8',newline='') as f:csv.writer(f,delimiter='\t',lineterminator='\n').writerow(t['columns'])
 catalog.append(dict(file=name+'.tsv',row_grain=t['grain'],data_scope=t['scope'],n_columns=len(t['columns'])))
with (OUT/'table_catalog.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(catalog[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(catalog)
config=dict(protocol='camp_discovery_to_cell_source_v1',status='TEMPLATE_NOT_RUN',cancer='REPLACE_ASSIGNED_CANCER',cohort='REPLACE_COHORT',run_id='REPLACE_UTC_RUN_ID',input_files=[],sample_design=dict(pair_evidence=None,unit=None,duplicate_policy=None,normal_tissue_definition=None),data_scale=dict(metabolite=None,RNA=None,availability_mask_semantics=None),min_pairs=8,min_correlation_n=8,metabolite_test=dict(profile='new_cohort_default_v1',test='two_sided_signed_rank',zero_method='wilcox',ties='average_ranks',all_zero_p=1,exact_signflip_max_nonzero=16,MC_signflip_max_nonzero=29,MC_B=99999,large_n_method='normal_tie_corrected_with_continuity',rounding_policy='predeclare_from_measurement_precision'),correlation=dict(method='spearman',p_method='two_sided_patient_permutation',B=9999,plus_one=True),bootstrap=dict(B=4000,unit='whole_patient_or_pair',ci='percentile95_pointwise',min_valid_fraction=.9),RNA=dict(input_type='MUST_SET',continuous_model='paired_t_on_documented_author_scale',count_model='separate_edgeR_paired_design_required'),multiplicity=dict(method='BH',family_denominator='all_evaluable_in_prespecified_family',missing_P_q='NA',families=['METAB_PAIRED_PRIMARY','METAB_PAIRED_AVAILABLE','REL_TUMOR_PRIMARY','REL_TUMOR_AVAILABLE','RNA_PAIRED_CURRENT']),workpool=dict(primary_metabolite_P_lt=.05,keep_all_q=True,no_RNA_or_correlation_gate_for_sc=True),single_cell=dict(min_cells_per_donor_type=20,min_labels_per_type=3,low_detection_flag=.01,bootstrap_B=1000,stability_tag=.8,rank_tie_rtol=1e-10,rank_tie_atol=1e-12,new_P_q=False,stop_after='expression_source',allow_new_umap_only_if_explicitly_scoped=True),master_seed=20260922,software_versions={},source_commit=None,planned_counts={},historical_results_immutable=True)
(OUT/'analysis_spec.template.json').write_text(json.dumps(config,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
(OUT/'README_CN.md').write_text('''# 空表与参数模板

配套规范：docs/CAMP_DISCOVERY_TO_CELL_SOURCE_SOP_CN.md。
这些文件只有表头，不是已完成分析或虚构结果；analysis_spec状态为TEMPLATE_NOT_RUN。复制到新癌种运行目录后填写已核实设计并冻结。不能直接把模板交付为DONE。

table_catalog.tsv说明每表一行代表什么、允许保存范围。SERVER_PRIVATE表填入真实ID/测量后只存server165；空模板可以公开。PUBLIC_AGGREGATE为群组统计。主统计前缀继承仓库templates/statistical_result.tsv。

效应主字段effect与effect_type配套；ci_lower/upper对应ci_estimand或明确的effect，不混用均值CI与秩相关CI。P/q不可测为NA并给reason；来源脚本/版本必须可追溯。相同统计量的补充列应核对一致。

临床字段、状态和计数按SOP填写，不能照抄BRCA计数。规范中的软件选择分支必须在analysis_spec中确定；特别是原始RNA计数不可直接套配对t。相同数据旧统计只复用，新q单独版本。

生成/校验空模板：python code/build_camp_source_sop_templates_v1.py。此命令不执行患者统计、不访问服务器、不上传真实数据。
''',encoding='utf-8')
for name,t in tables.items():
 rows=list(csv.reader((OUT/(name+'.tsv')).open(encoding='utf-8'),delimiter='\t'));assert rows==[t['columns']]
val=dict(status='TEMPLATES_VALIDATED_NOT_ANALYSIS',tables=len(tables),empty_tables_only=True,unique_headers=True,statistical_prefix_matches_repository=True,config_json_parseable=True,patient_data_included=False)
(OUT/'validation.json').write_text(json.dumps(val,indent=2)+'\n')
(OUT/'.gitattributes').write_text('* -text\n',encoding='utf-8')
paths=[p for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='checksums.tsv']
with (OUT/'checksums.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['file','sha256']);w.writerows((p.name,hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths)
print(json.dumps(val))
