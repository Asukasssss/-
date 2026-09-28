"""Freeze public aggregate inputs for the presentation. No patient data or new tests."""
from pathlib import Path
import subprocess,hashlib,json
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
RUN='20260928T150000Z_ppt_redraw_v1'
OUT=ROOT/'results/BRCA/07_INTEGRATION'/RUN
SRC=OUT/'sources';SRC.mkdir(parents=True,exist_ok=True)
BR='4b2a02b5cdeb999ec94f4558f93ce2f86a743a54'
CO='e18818697125f457fbd1c28db7c4dbb90f3de72f'
PD='9239959a17175222778301fdea6a094f9509242e'
PR='8fe4d512390c5005da36be4920e20ae2517d214a'
GB='2174aa965e29cc0fbb10925bc2c217198960829e'
CC='02593f745fc966c975ab1a25e3e7ce74f770efc6'
ST='a395909c8f4a6234cd8b2faec9db344ceb4aa056'
spec={
'brca_met':(BR,'results/BRCA/07_INTEGRATION/20260922T140000Z_report_v1/figure_data/02_metabolites.tsv'),
'coad_met':(CO,'results/COAD/07_INTEGRATION/20260922T150804Z_source_report_v2/figure_data/02_metabolites.tsv'),
'pdac_met':(PD,'results/PDAC/01_CAMP/20260922T121000Z_sop_v3_discovery_complete/metabolite_paired.tsv'),
'prad_met':(PR,'results/PRAD/07_INTEGRATION/20260922T142000Z_discovery_v1/metabolite_all361_integrated.tsv'),
'gbm_met':(GB,'results/GBM/07_INTEGRATION/20260925T145000Z_integration_v1/all_metabolites_current_view.tsv'),
'ccrcc_met':(CC,'results/ccRCC/07_INTEGRATION/20260925T163000Z_integration_v1/figures/02_metabolite_discovery_source.tsv'),
'brca_rel':(BR,'results/BRCA/07_INTEGRATION/20260921T154000Z_mapping262_integration_v1/relations262_comparison.tsv'),
'brca_gene':(BR,'results/BRCA/07_INTEGRATION/20260921T154000Z_mapping262_integration_v1/genes156_all_history_comparison.tsv'),
'coad_rel':(CO,'results/COAD/07_INTEGRATION/20260922T151000Z_source_contract_integration_v2/candidate_relations.tsv'),
'coad_gene':(CO,'results/COAD/07_INTEGRATION/20260922T151000Z_source_contract_integration_v2/candidate_genes.tsv'),
'pdac_rel':(PD,'results/PDAC/07_INTEGRATION/20260922T131900Z_author_identity_v4/candidate_relations_integrated.tsv'),
'pdac_gene':(PD,'results/PDAC/07_INTEGRATION/20260922T131900Z_author_identity_v4/candidate_genes_integrated.tsv'),
'prad_rel':(PR,'results/PRAD/07_INTEGRATION/20260922T142000Z_discovery_v1/candidate_relations_integrated.tsv'),
'prad_gene':(PR,'results/PRAD/07_INTEGRATION/20260922T142000Z_discovery_v1/candidate_genes_integrated.tsv'),
'brca_sc':(BR,'results/BRCA/06_EXTERNAL/20260923T040012Z_epithelial_paired20_v1/results.tsv'),
'coad_sc':(CO,'results/COAD/06_EXTERNAL/20260925T142407Z_lypla1_cell_pooled_v1/results.tsv'),
'pdac_sc':(PD,'results/PDAC/06_EXTERNAL/20260928T031629Z_gse202051_lypla1_v1/results.tsv'),
'pdac_expr':(PD,'results/PDAC/06_EXTERNAL/20260928T031629Z_gse202051_lypla1_v1/expression_summary.tsv'),
'prad_expr':(PR,'results/PRAD/06_EXTERNAL/20260928T031839Z_lypla1_umap_v2/epithelial_summary.tsv'),
'prad_positive':(PR,'results/PRAD/06_EXTERNAL/20260928T090606Z_lypla1_positive_summary_v1/positive_cell_summary.tsv'),
'brca_spatial':(BR,'results/BRCA/06_EXTERNAL/20260928T060000Z_cta_lypla1_v1/aggregate_effects.tsv'),
'pdac_spatial':(PD,'results/PDAC/06_EXTERNAL/20260928T094758Z_geomx_bell2025_lypla1_v1/results.tsv'),
'prad_spatial':(PR,'results/PRAD/06_EXTERNAL/20260928T101700Z_lypla1_erickson_v1/contrast_summary.tsv'),
'prad_regions':(PR,'results/PRAD/06_EXTERNAL/20260928T101700Z_lypla1_erickson_v1/region_summary.tsv'),
'brca_states':(ST,'results/BRCA/06_EXTERNAL/20260928T131206Z_lypla1_cellstates_v1/state_summary.tsv'),
'brca_companions':(ST,'results/BRCA/06_EXTERNAL/20260928T131206Z_lypla1_cellstates_v1/cross_cohort_companions.tsv'),
'brca_delta':(BR,'results/BRCA/03_PATIENT/20260925T091556Z_lypla1_paired_delta_v1/results.tsv'),
}
manifest=[]
for name,(sha,path) in spec.items():
 b=subprocess.check_output(['git','show',sha+':'+path],cwd=ROOT)
 (SRC/(name+'.tsv')).write_bytes(b)
 manifest.append(dict(id=name,git_commit=sha,git_path=path,sha256=hashlib.sha256(b).hexdigest(),bytes=len(b),url=f'https://github.com/Asukasssss/-/blob/{sha}/{path}'))
pd.DataFrame(manifest).to_csv(OUT/'source_manifest.tsv',sep='\t',index=False)
(OUT/'analysis_spec.json').write_text(json.dumps(dict(version='ppt_redraw_v1',scope='Public aggregate reuse; four-cancer LYPLA1 focus; six-cancer CAMP background',new_statistics=False,new_embedding=False,source_matrices='server165 only',nominal_threshold=.05,source_manifest='source_manifest.tsv'),indent=2),encoding='utf8')
print('Frozen',len(manifest),'aggregate source tables')
for name in ['brca_met','coad_met','pdac_met','gbm_met','ccrcc_met']:
 d=pd.read_csv(SRC/(name+'.tsv'),sep='\t'); print(name,len(d),'types',d.get('analysis_type',pd.Series()).unique(),'cohorts',d.get('cohort',pd.Series()).unique())
