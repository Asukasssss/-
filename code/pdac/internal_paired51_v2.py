"""Server-only PDAC 21-author-unit associations and 11-pair RNA background."""
import argparse,hashlib,json,platform,traceback
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from patient_first_round import compute,leave_one_out,fixed_bh,validate_rna_labels
from paired_metabolites_v1 import pair_join,paired_stat
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
SRC=ROOT/'data/candidates/camp_primary_tissue_multicancer'
MP=SRC/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv'
XP=SRC/'processed_metabolomics/PreprocessedData_PDAC.xlsx'
RP=SRC/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed/GSE62452.hugene10st.gene_symbol.csv'
AP=ROOT/'results/collaborative/PDAC/B/20260921T120900Z_paired_source_v1/GSE62452_paired_sample_information.xlsx'
EXPECTED={str(MP):'7a5332cfdeba3cc055d79ef223c95a3c8f18ecf6a7058fafd72b390aec62dce7',str(XP):'e5f408848824394d6be875d7284b77d0a0aabd433d0f2db8ebe01886c29a7ecd',str(RP):'d2fb75b102cb5adf7264b66e6664e5b2b0a07241b3f9901f125c78959da2e8ba',str(AP):'0f16963287abb1bf7bf73dd5a2738f3bdccfb0344b2db3db72598199432ac866'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def validate_sources(expected,required):
 if not isinstance(expected,dict) or not required or set(expected)!=set(map(str,required)):raise ValueError('Exact actual input-path manifest required')
 for p,h in expected.items():
  if not isinstance(h,str) or len(h)!=64 or sha(p)!=h:raise ValueError('Source hash mismatch: '+p)
def resolve_genes(genes,labels,aliases):
 labels=set(labels);out={};reasons={}
 for g in genes:
  found=[g] if g in labels else sorted(set(aliases.get(g,[]))&labels)
  out[g]=found[0] if len(found)==1 else None
  reasons[g]='EXACT_CANONICAL' if found==[g] else 'UNIQUE_REVIEWED_ALIAS' if len(found)==1 else 'MISSING_OR_AMBIGUOUS'
 counts={s:list(out.values()).count(s) for s in set(out.values()) if s}
 for g,s in list(out.items()):
  if s and counts[s]>1:out[g]=None;reasons[g]='LABEL_COLLISION_BETWEEN_TARGET_GENES'
 return out,reasons
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--code-commit',required=True);args=ap.parse_args()
 run=Path(__file__).resolve().parent;assert run.parent==ROOT/'results/collaborative/PDAC/B'
 with (run/'.running').open('x') as f:f.write('PDAC internal paired51 v2')
 try:
  out=run/'public';out.mkdir(exist_ok=False)
  expected=json.loads((run/'input_hashes.json').read_text());validate_sources(expected,[MP,XP,RP,AP]);assert expected==EXPECTED
  rel=pd.read_csv(run/'planned_relations.tsv',sep='\t',keep_default_na=False);panel=json.loads((run/'panel.json').read_text());aliases=json.loads((run/'gene_aliases.json').read_text());genes=panel['genes']
  assert len(rel)==715 and rel.relation_id.is_unique and rel.groupby('mapping_tier').size().to_dict()=={'CONDITIONAL':358,'DIRECT':357}
  assert len(genes)==538 and set(genes)==set(rel.gene)
  joined,pairs=pair_join(pd.read_csv(MP,dtype=str),pd.read_excel(AP,dtype=str));t=joined[(joined.TN=='Tumor')&joined['Paired Sample ID'].notna()].sort_values('Paired Sample ID').copy()
  assert len(t)==21 and t['Paired Sample ID'].is_unique and len(pairs)==11
  paired_ids=set(pairs.author_pair_id);tp=joined[(joined.TN=='Tumor')&joined['Paired Sample ID'].isin(paired_ids)].sort_values('Paired Sample ID');npair=joined[(joined.TN=='Normal')&joined['Paired Sample ID'].isin(paired_ids)].sort_values('Paired Sample ID')
  assert len(tp)==len(npair)==11 and tp['Paired Sample ID'].tolist()==npair['Paired Sample ID'].tolist()
  joined.to_csv(run/'private_author_join.tsv',sep='\t',index=False)
  rna=pd.read_csv(RP,index_col=0);validate_rna_labels(rna)
  met=pd.read_excel(XP,sheet_name='metabo_imputed_filtered_Tumor',index_col=0);raw=pd.read_excel(XP,sheet_name='data',index_col=0)
  for df in [met,raw]:validate_rna_labels(df)
  assert joined.RNAID.isin(rna.columns).all() and t.MetabID.isin(met.columns).all() and t.MetabID.isin(raw.columns).all()
  assert set(rel.feature_name)<=set(met.index)&set(raw.index)
  resolved,reasons=resolve_genes(genes,rna.index,aliases)
  pd.DataFrame([{'gene':g,'rna_label':resolved[g],'reason':reasons[g]} for g in genes]).to_csv(out/'gene_resolution.tsv',sep='\t',index=False,na_rep='NA')
  rr=[];cache={}
  for i,r in rel.iterrows():
   x=met.loc[r.feature_name,t.MetabID].to_numpy(float);avail=np.isfinite(raw.loc[r.feature_name,t.MetabID].to_numpy(float));label=resolved[r.gene]
   y=rna.loc[label,t.RNAID].to_numpy(float) if label else np.full(len(t),np.nan)
   for mode in ['primary','availability']:
    mask=np.isfinite(x)&np.isfinite(y)
    if mode=='availability':mask&=avail
    a=x[mask];b=y[mask];key=hashlib.sha256(r.relation_id.encode()+mask.tobytes()+a.tobytes()+b.tobytes()+b'PDAC_internal_v2').hexdigest();seed=int(key[:8],16);reused=key in cache
    if not reused:
     z=compute(a,b,seed);z['reason']='Author-confirmed distinct pair-key tumor units;exploratory association' if z['status']=='DONE' else reasons[r.gene] if not label else z['reason'];cache[key]=(z,leave_one_out(a,b))
    z,loo=cache[key]
    rr.append({'cancer':'PDAC','cohort':'CAMP_PDAC_GSE62452','stage_id':'03_PATIENT','run_id':run.name,'analysis_version':'PDAC_internal_paired51_v2','analysis_type':mode,'metabolite_key':r.metabolite_key,'metabolite_name':r.feature_name,'gene':r.gene,'unit':'author_distinct_pair_key_tumor','n':len(a),'n_reference':21,'effect_type':'Spearman_rho','effect':z.get('rho',np.nan),'ci_lower':z.get('lo',np.nan),'ci_upper':z.get('hi',np.nan),'p_value':z.get('p',np.nan),'q_value':np.nan,'test_family':r.mapping_tier+'_'+mode,'family_n_evaluable':0,'status':z['status'],'reason':z['reason'],'source_id':'CAMP_plus_GEO_author_pair_table','relation_id':r.relation_id,'mapping_tier':r.mapping_tier,'rna_label':label,'n_preimputation_available':int(avail.sum()),'available_fraction':float(avail.mean()),'same_input_reused':reused,'input_hash':key,'seed':seed,'bootstrap_valid':z.get('bootstrap_valid',0),**loo})
   if i%100==0:print('relationships',i+1,'/',len(rel),flush=True)
  results=pd.DataFrame(rr)
  for family,d in results.groupby('test_family'):
   p=[None if pd.isna(v) else v for v in d.p_value];results.loc[d.index,'q_value']=fixed_bh(p);results.loc[d.index,'family_n_evaluable']=d.p_value.notna().sum();results.loc[d.index,'family_n_planned']=len(d)
  results.to_csv(out/'associations.tsv',sep='\t',index=False,na_rep='NA')
  er=[]
  for g in genes:
   lab=resolved[g];a=rna.loc[lab,tp.RNAID].to_numpy(float) if lab else np.full(11,np.nan);b=rna.loc[lab,npair.RNAID].to_numpy(float) if lab else np.full(11,np.nan);mask=np.isfinite(a)&np.isfinite(b);delta=(a-b)[mask];seed=int(hashlib.sha256(('PDAC_RNA_paired51_v2|'+g).encode()).hexdigest()[:8],16);z=paired_stat(delta,seed)
   er.append({'cancer':'PDAC','cohort':'CAMP_PDAC_GSE62452','stage_id':'03_PATIENT','run_id':run.name,'analysis_version':'PDAC_RNA_paired51_v2','analysis_type':'paired_RNA','metabolite_key':'NA','metabolite_name':'NA','gene':g,'unit':'author_confirmed_patient_pair','n':len(delta),'n_reference':11,'effect_type':'paired_mean_tumor_minus_normal_author_RNA_scale','effect':z.get('mean_delta',np.nan),'ci_lower':z.get('ci_lower',np.nan),'ci_upper':z.get('ci_upper',np.nan),'p_value':z.get('p_value',np.nan),'q_value':np.nan,'test_family':'all538_candidate_genes_paired_RNA','family_n_evaluable':0,'status':z['status'],'reason':z['reason'],'source_id':'CAMP_plus_GEO_author_pair_table','rna_label':lab,'median_delta':z.get('median_delta',np.nan),'rank_biserial':z.get('rank_biserial',np.nan),'n_up':int((delta>0).sum()),'n_down':int((delta<0).sum()),'n_equal':int((delta==0).sum()),'seed':seed,'family_n_planned':len(genes)})
  expression=pd.DataFrame(er);expression.q_value=fixed_bh([None if pd.isna(x) else x for x in expression.p_value]);expression.family_n_evaluable=int(expression.p_value.notna().sum());expression.to_csv(out/'paired_RNA.tsv',sep='\t',index=False,na_rep='NA')
  pd.DataFrame([{'scope':s,'count':n,'status':'DONE','basis':b} for s,n,b in [('CAMP_tumor_specimens',27,'Exact CAMP PDAC mapping'),('CAMP_normal_specimens',12,'Exact CAMP PDAC mapping'),('tumor_author_distinct_pair_key',21,'GSM join to explicit author pairing workbook;all unique'),('tumor_missing_author_pair_key_excluded',6,'Not assumed independent'),('complete_CAMP_metabolite_and_RNA_pairs',11,'Explicit author T/N pair groups'),('tissue_label_conflicts',0,'All39 checked;observed zero'),('unmatched_CAMP_RNA_columns',0,'All39 exact matrix column matches')]]).to_csv(out/'sample_identity_audit.tsv',sep='\t',index=False)
  summary={'tumor_units':21,'RNA_pairs':11,'relations':715,'genes':538,'paired_RNA_evaluable':int(expression.p_value.notna().sum()),'paired_RNA_P_lt005':int((expression.p_value<.05).sum()),'paired_RNA_q_lt005':int((expression.q_value<.05).sum()),'families':{k:{'planned':len(d),'evaluable':int(d.p_value.notna().sum()),'P_lt005':int((d.p_value<.05).sum()),'q_lt005':int((d.q_value<.05).sum())} for k,d in results.groupby('test_family')}}
  save(out/'summary.json',summary);save(out/'analysis_spec.json',{'code_commit':args.code_commit,'protocol':'docs/PDAC/PAIRED51_INTERNAL_V2_LOCK.md','planned_panel':panel,'association':{'minimum_n':8,'permutations':9999,'bootstrap':4000,'fixed_BH_families':'DIRECT and CONDITIONAL separately,primary and availability separately'},'paired_RNA':{'minimum_pairs':6,'exact_signed_rank':True,'bootstrap':10000,'fixed_BH':538},'no_new_imputation_or_transform':True,'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__}})
  sources=[MP,XP,RP,AP,*[run/n for n in ['planned_relations.tsv','panel.json','gene_aliases.json','input_hashes.json','patient_first_round.py','paired_metabolites_v1.py']],Path(__file__)]
  pd.DataFrame([{'path':str(p),'sha256':sha(p)} for p in sources]).to_csv(out/'source_manifest.tsv',sep='\t',index=False)
  assert len(results)==1430 and len(expression)==538 and not results.duplicated(['relation_id','analysis_type']).any()
  save(out/'validation.json',{'status':'PASS','exact_four_input_hash_gate':True,'RNA_labels_complete_unique':True,'author_pair_key_unique':21,'paired_RNA_patient_alignment':11,'all_locked_rows_retained':True,'patient_level_rows_exported':False,'not_verified':['Clinical covariates','Original frozen CAMP contrast code','Chemical reidentification','Causal mechanism','External matched-omics replication']})
  save(run/'DONE.json',summary);(run/'.running').unlink();print(json.dumps(summary),flush=True)
 except Exception:(run/'FAILED.txt').write_text(traceback.format_exc());raise
if __name__=='__main__':main()
