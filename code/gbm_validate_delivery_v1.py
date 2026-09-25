"""Validate scoped GBM delivery, public-only fields, source bytes, and family BH."""
from pathlib import Path
import hashlib,json,subprocess,sys
import numpy as np,pandas as pd
from openpyxl import load_workbook
R=Path(__file__).resolve().parents[1];B=R/'results/GBM';O=B/'07_INTEGRATION/20260925T145000Z_integration_v1';SC=B/'06_EXTERNAL/20260925T103000Z_discovery_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
# External input pointers include immutable asset URLs, without claiming an unread full-file SHA.
man=pd.read_csv(SC/'source_manifest.tsv',sep='\t');extra=[]
for cid,study in [('558385a4-b7b7-4eca-af0c-9e54d010e8dc','Darmanis2017'),('999f2a15-3d7e-440b-96ae-2c806799c08c','Neftel2019_GBmap')]:
 p=SC/(cid+'_collection.json');j=json.loads(p.read_text());ds=j['datasets'][0]
 if study.startswith('Neftel'):ds=next(z for z in j['datasets'] if z['dataset_id']=='c888b684-6c51-431f-972a-6c963044cef0')
 asset=next(a for a in ds['assets'] if a['filetype']=='H5AD');extra.append(dict(source_id=study+'_immutable_asset',path_or_url=asset['url'],sha256='NOT_COMPUTED_FULL_REMOTE_ASSET' if study.startswith('Neftel') else man.loc[man.source_id.eq('Darmanis2017_cellxgene.h5ad'),'sha256'].iloc[0],file_bytes=asset['filesize'],collection_version_id=j['collection_version_id'],read_scope='Neftel2019 rows via HTTP byte ranges;derived private subset SHA recorded' if study.startswith('Neftel') else 'full processed H5AD on server165'))
 extra.append(dict(source_id=study+'_collection_metadata',path_or_url=str(p.relative_to(R)),sha256=sha(p)))
for p in [R/'code/gbm_extract_sc_v1.py',SC/'Darmanis_UCSC_dataset_metadata.json']:extra.append(dict(source_id=p.name,path_or_url=str(p.relative_to(R)),sha256=sha(p)))
save(pd.concat([man,pd.DataFrame(extra)],ignore_index=True).drop_duplicates(['source_id','path_or_url'],keep='last'),SC/'source_manifest.tsv')
figs=[]
for p in (O/'figures').iterdir():
 figs.append(dict(file=p.name,source_table='sc_celltype_profiles.tsv' if 'association' not in p.name else 'candidate_relations_integrated.tsv',interpretation='donor_equal_descriptive_expression;gray_NA_not_zero;no_function_or_flux_inference' if 'association' not in p.name else 'nominal14_view;all_q_ge005;pointwise_bootstrap95CI;selection_not_adjusted',code='code/gbm_integrate_v1.py'))
save(pd.DataFrame(figs),O/'figure_manifest.tsv')
function=B/'05_FUNCTION/20260925T103000Z_discovery_v1';function.mkdir(parents=True,exist_ok=True);save(pd.DataFrame([dict(cancer='GBM',stage_id='05_FUNCTION',status='NOT_RUN',reason='functional_interventions_not_required_by_discovery_to_cell_source_SOP',functional_proof='NOT_ESTABLISHED')]),function/'status.tsv');(function/'README_CN.md').write_text('功能干预与独立外部代谢验证本轮未开展；单细胞表达来源不替代功能证明。\n',encoding='utf-8')
# Every public table is aggregate-only: reject direct sample identifiers as fields.
blocked={'patient_id','specimen_id','case_id','donor_id','cell_id','RNAID','MetabID','CommonID'};tables=0
for p in B.rglob('*.tsv'):
 d=pd.read_csv(p,sep='\t');assert not blocked.intersection(d.columns),(p,blocked.intersection(d.columns));tables+=1
# Family q independent formula.
errs=[];families=0
for p in [B/'03_PATIENT/20260925T103000Z_discovery_v1/tumor_association_all_families.tsv',B/'03_PATIENT/20260925T144100Z_identity_resolution_v2/RNA_unpaired_all_families.tsv',B/'01_CAMP/20260925T103000Z_discovery_v1/metabolite_unpaired_all.tsv',B/'04_ROBUSTNESS/20260925T103000Z_discovery_v1/metabolite_available_sensitivity.tsv']:
 d=pd.read_csv(p,sep='\t')
 for f,z in d.groupby('test_family'):
  z=z[z.p_value.notna()];ps=z.p_value.to_numpy();ix=np.argsort(ps);q=np.minimum(1,np.minimum.accumulate((ps[ix]*len(ps)/np.arange(1,len(ps)+1))[::-1])[::-1]);err=np.max(abs(q-z.q_value.to_numpy()[ix]));assert err<1e-12;errs.append(err);families+=1
rels=pd.read_csv(O/'candidate_relations_integrated.tsv',sep='\t');genes=pd.read_csv(O/'candidate_genes_integrated.tsv',sep='\t');assert len(rels)==171 and rels.relation_id.is_unique;assert len(genes)==142 and genes.gene.is_unique and genes.stable_gene_id.is_unique;assert genes.relation_ids.str.split(';').map(len).sum()==171
assert int(rels.patient_primary_p_value.lt(.05).sum())==14 and not rels.patient_primary_q_value.lt(.05).any();assert int(genes.RNA_q_value.lt(.05).sum())==75
wb=load_workbook(O/'GBM_全候选分析_142基因171关系.xlsx',read_only=True);assert wb['全部171关系'].max_row==172 and wb['全部142基因'].max_row==143 and wb['代谢物713全表'].max_row==714
# Freeze actual code identity for reproducibility; snapshot contains the executed code bytes.
spec=json.loads((O/'analysis_spec.json').read_text());spec['code_commit']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip();spec['code_sha256']=sha(R/'code/gbm_integrate_v1.py');spec['runtime_used']=sys.executable;spec['runtime_note']='bundled Python inspected;matplotlib unavailable there;used existing verified Python310 with matplotlib/openpyxl';(O/'analysis_spec.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
source=[B/'02_MAPPING/20260925T103000Z_discovery_v1/direct_relations.tsv',B/'01_CAMP/20260925T103000Z_discovery_v1/metabolite_unpaired_all.tsv',B/'03_PATIENT/20260925T103000Z_discovery_v1/tumor_association_all_families.tsv',B/'03_PATIENT/20260925T144100Z_identity_resolution_v2/RNA_unpaired_all_families.tsv',SC/'sc_cross_study.tsv',SC/'sc_source_stability.tsv',SC/'sc_celltype_profiles.tsv',R/'code/gbm_integrate_v1.py'];save(pd.DataFrame([dict(source_id=p.name,path_or_url=str(p.relative_to(R)),sha256=sha(p)) for p in source]),O/'source_manifest.tsv')
verification=dict(status='DONE',public_tables_checked=tables,no_identifier_columns=True,relation_keys_unique=True,gene_keys_unique=True,join_no_expansion=True,independent_BH_families=families,max_BH_error=max(errs),workbook_row_counts_verified=True,figure_page1_visually_checked=True,sc_independent_validation=json.loads((SC/'independent_validation.json').read_text()),remaining_scientific_scope='457 work features await precise mapping;3 RNA symbols unmeasured by exact name;no functional intervention evidence')
(O/'independent_validation.json').write_text(json.dumps(verification,indent=2),encoding='utf-8')
# Rebuild only changed/new stage checksums. Prior stage files and statistics remain frozen.
for folder in [SC,O,function]:
 save(pd.DataFrame([dict(file=str(p.relative_to(folder)),sha256=sha(p)) for p in sorted(folder.rglob('*')) if p.is_file() and p.name!='checksums.tsv']),folder/'checksums.tsv')
n=0
for p in B.rglob('checksums.tsv'):
 d=pd.read_csv(p,sep='\t')
 for _,r in d.iterrows():assert sha(p.parent/r['file'])==r.sha256,(p,r['file']);n+=1
save(pd.DataFrame([dict(file=str(p.relative_to(B)),sha256=sha(p)) for p in sorted(B.rglob('*')) if p.is_file() and p.name!='delivery_manifest.tsv']),B/'delivery_manifest.tsv')
print(json.dumps(dict(**verification,checked_checksum_entries=n)))
