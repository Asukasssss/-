"""Server-only follow-up; export full aggregate catalogs, never patient rows."""
import argparse,json,platform,traceback
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from patient_statistics import bh,stable_seed,association
from covariate_statistics import bootstrap_plan,adjusted
from run_patient_v1 import read,table,dump,sha

def main():
 p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,required=True);p.add_argument('--run-dir',type=Path,required=True);p.add_argument('--code-commit',required=True);a=p.parse_args()
 root,o=a.project_root.resolve(),a.run_dir.resolve();base=root/'results/collaborative/COAD/B'
 assert o.parent==base and (o/'.running').read_text().strip()==o.name and not (o/'analysis_spec.json').exists()
 res=base/'20260919T115842Z_identity_units_v1';old=base/'20260919T122034Z_patient_v1'
 spec=json.loads((o/'locked_spec.json').read_text());assert json.loads((o/'kernel_validation.json').read_text())['status']=='PASS'
 clinical=res/'private_clinical_join.tsv';assert sha(clinical)==spec['clinical_sha256']
 c=pd.read_csv(clinical,sep='\t',dtype=str,keep_default_na=False);c=c[c.TN.eq('Tumor') & c.stage.isin(spec['categories']['stage'])].copy()
 assert len(c)==33 and c.individual.is_unique
 data=root/'data/candidates/camp_primary_tissue_multicancer';mp=data/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv'
 mapping=pd.read_csv(mp,dtype=str);mapping=mapping[mapping.Dataset.eq('COAD') & mapping.TN.eq('Tumor')]
 rp=data/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/mapping.RNAFile.iloc[0]
 xp=data/'processed_metabolomics'/mapping.MetabFile.iloc[0]
 rna=pd.read_csv(rp,index_col=0);met=pd.read_excel(xp,sheet_name=mapping.MetabFile_sheet.iloc[0],index_col=0)
 for matrix in [rna,met]:
  matrix.index=matrix.index.astype(str);matrix.columns=matrix.columns.astype(str);assert matrix.index.is_unique and matrix.columns.is_unique
 candidates=read(res/'candidate_eligibility.tsv');assert len(candidates)==974 and sum(r['in_prespecified_family']=='True' for r in candidates)==674
 prior=read(old/'results.tsv');key=lambda r:(r['metabolite_key'],r['gene']);oldrows={key(r):r for r in prior};assert len(oldrows)==974
 priorvalidation=json.loads((old/'validation.json').read_text());assert sha(old/'results.tsv')==priorvalidation['public_output_sha256']['results.tsv']
 priorinput=json.loads((old/'analysis_spec.json').read_text())['input_hashes']
 for f in [clinical,rp,xp,res/'candidate_eligibility.tsv']:assert sha(f)==priorinput[str(f.relative_to(root))]
 inputs=[clinical,rp,xp,mp,res/'candidate_eligibility.tsv',old/'results.tsv',o/'locked_spec.json',o/'statistical_result_template.tsv']+sorted(o.glob('*.py'))
 hashes={str(f):sha(f) for f in inputs}
 dump(o/'analysis_spec.json',dict(spec,run_id=o.name,code_commit=a.code_commit,input_hashes=hashes,software={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__}))
 prefix=(o/'statistical_result_template.tsv').read_text().strip().split('\t')
 extra=['human_gene_id','rna_symbol','original_effect','original_q','family_n_planned','model','nuisance_columns','design_rank','residual_df','min_category_n','n_excluded','bootstrap_valid','bootstrap_invalid','bootstrap_seed','ci_status','permutations','random_seed','ols_max_abs_error','loo_valid','loo_invalid','loo_min','loo_max','loo_max_abs_delta','loo_any_sign_change','m0_effect','m0_p','m0_q','m0_n','m0_cc_effect','m0_cc_p','m0_cc_q','m0_cc_n','m0_cc_origin','delta_adjustment','delta_case_selection']
 fields=prefix+extra;all_results=[];all_cc=[];summaries={};designs={};bootcache={}
 age=pd.to_numeric(c.age,errors='coerce').to_numpy(dtype=float)
 for model,catnames in spec['models'].items():
  mask=np.isfinite(age)
  for k in catnames:mask &=c[k].isin(spec['categories'][k]).to_numpy()
  z=np.column_stack([np.ones(len(c)),age-69]+[(c[k].eq(lev)).to_numpy(dtype=float) for k in catnames for lev in spec['categories'][k][1:]])
  designs[model]=(mask,z,catnames)
 for model,(covmask,zall,catnames) in designs.items():
  rows=[];ccrows=[]
  for i,r in enumerate(candidates):
   oldr=oldrows[key(r)]
   row={k:None for k in fields};row.update({k:oldr.get(k) for k in ['cancer','cohort','metabolite_key','metabolite_name','gene','human_gene_id','rna_symbol','unit','source_id','original_effect','original_q']})
   row.update(stage_id='04_ROBUSTNESS',run_id=o.name,analysis_version=spec['version'],analysis_type=model,model=model,effect_type='partial_rank_rho',test_family=spec['version']+'_'+model,family_n_planned=674,m0_effect=oldr['effect'],m0_p=oldr['p_value'],m0_q=oldr['q_value'],m0_n=oldr['n'])
   cc=dict(row,analysis_type=model+'_M0_CC',effect_type='Spearman_rho',test_family=spec['version']+'_'+model+'_M0_CC')
   if r['in_prespecified_family']!='True':
    row.update(status='NEEDS_REVIEW',reason='OUTSIDE_LOCKED_FAMILY',test_family='NOT_IN_CURRENT_FAMILY',family_n_planned=None);cc.update(status='NEEDS_REVIEW',reason='OUTSIDE_LOCKED_FAMILY',test_family='NOT_IN_CURRENT_FAMILY',family_n_planned=None)
   elif r['rna_symbol']=='NA':
    row.update(status='NOT_EVALUABLE',reason='RNA_SYMBOL_UNAVAILABLE');cc.update(status='NOT_EVALUABLE',reason='RNA_SYMBOL_UNAVAILABLE')
   else:
    x=pd.to_numeric(met.loc[r['metabolite_name'],c.MetabID],errors='coerce').to_numpy(dtype=float)
    y=pd.to_numeric(rna.loc[r['rna_symbol'],c.RNAID],errors='coerce').to_numpy(dtype=float)
    mask=covmask & np.isfinite(x) & np.isfinite(y);xx=x[mask];yy=y[mask];z=zall[mask];n=int(mask.sum());rank=int(np.linalg.matrix_rank(z)) if n else 0
    minn=min(int(c.loc[mask,k].eq(lev).sum()) for k in catnames for lev in spec['categories'][k])
    row.update(n=n,n_excluded=33-n,nuisance_columns=z.shape[1],design_rank=rank,residual_df=n-rank-1,min_category_n=minn)
    cc.update(n=n,n_excluded=33-n)
    if mask.all() and oldr['status']=='DONE':
     assert abs(stats.spearmanr(xx,yy).statistic-float(oldr['effect']))<1e-12
     cc.update({k:oldr[k] for k in ['n','effect','ci_lower','ci_upper','p_value','q_value','bootstrap_valid','permutations','random_seed']});cc.update(status='DONE',reason='IDENTICAL_INPUT_PRIMARY_REUSED',m0_cc_origin='IDENTICAL_INPUT_PRIMARY_REUSED')
    else:
     cc.update(association(xx,yy,stable_seed(spec['version']+'|'+model+'|M0_CC|'+r['metabolite_key']+'|'+r['human_gene_id']),minimum_n=8));cc['m0_cc_origin']='NEW_SUBSET_CALCULATION'
    reason=''
    if n<spec['minimum_n']:reason='N_BELOW_MINIMUM'
    elif rank<z.shape[1]:reason='NUISANCE_DESIGN_RANK_DEFICIENT'
    elif n-rank-1<spec['minimum_full_model_residual_df']:reason='INSUFFICIENT_RESIDUAL_DF'
    elif minn<spec['minimum_category_n']:reason='SPARSE_SOURCE_CATEGORY_LT3'
    if reason:row.update(status='NOT_EVALUABLE',reason=reason)
    else:
     maskid=''.join('1' if v else '0' for v in mask);cachekey=(model,maskid)
     bs=stable_seed(spec['version']+'|bootstrap|'+model+'|'+maskid)
     if cachekey not in bootcache:bootcache[cachekey]=bootstrap_plan(z,bs,spec['bootstrap_replicates'])
     seed=stable_seed('|'.join([spec['version'],model,r['metabolite_key'],r['human_gene_id']]))
     row.update(adjusted(xx,yy,z,seed,spec,bootcache[cachekey]));row['bootstrap_seed']=bs
    if cc['status']=='DONE':
     row.update(m0_cc_effect=cc['effect'],m0_cc_p=cc['p_value'],m0_cc_n=cc['n'],m0_cc_origin=cc['m0_cc_origin'])
     row['delta_case_selection']=float(cc['effect'])-float(oldr['effect']) if oldr['effect']!='NA' else None
     if row['status']=='DONE':row['delta_adjustment']=row['effect']-float(cc['effect'])
   rows.append(row);ccrows.append(cc)
   if (i+1)%100==0:print(model,i+1,'/974',flush=True)
  for rs in [rows,ccrows]:
   good=[r for r in rs if r['status']=='DONE'];q=bh([float(r['p_value']) for r in good])
   for r,v in zip(good,q):r['q_value']=float(v)
   for r in rs:
    if r['test_family']!='NOT_IN_CURRENT_FAMILY':r['family_n_evaluable']=len(good)
  for row,cc in zip(rows,ccrows):
   assert key(row)==key(cc)
   if cc['status']=='DONE':row['m0_cc_q']=cc['q_value']
  good=[r for r in rows if r['status']=='DONE'];summaries[model]={'planned':674,'evaluable':len(good),'p_lt_005':sum(r['p_value']<.05 for r in good),'q_lt_005':sum(r['q_value']<.05 for r in good),'min_q':min((r['q_value'] for r in good),default=None),'n_distribution':dict(Counter(str(r['n']) for r in good)),'unavailable_reasons':dict(Counter(r['reason'] for r in rows if r['status']=='NOT_EVALUABLE')),'m0_cc_reused':sum(r.get('m0_cc_origin')=='IDENTICAL_INPUT_PRIMARY_REUSED' for r in ccrows),'loo_sign_change':sum(r.get('loo_any_sign_change',False) for r in good)}
  order=lambda r:(r['metabolite_name'],r['metabolite_key'],r['gene'],r['human_gene_id'])
  table(o/(model+'.tsv'),fields,sorted(rows,key=order));table(o/(model+'_M0_CC.tsv'),fields,sorted(ccrows,key=order))
  all_results+=rows;all_cc+=ccrows
 for f,h in hashes.items():assert sha(Path(f))==h
 assert len(all_results)==3*974
 for r in all_results:
  oldr=oldrows[key(r)];assert r['original_q']==oldr['original_q'] and r['original_effect']==oldr['original_effect']
  if r['status']!='DONE':assert r['p_value'] is None and r['q_value'] is None
 summary={'run_id':o.name,'primary_n':33,'models':summaries,'new_patient_statistics':True,'purity':'NOT_RUN','technical_batch':'NOT_RUN','external_validation':'NOT_RUN'}
 dump(o/'summary.json',summary);table(o/'source_manifest.tsv',['source_path','sha256'],[{'source_path':f,'sha256':h} for f,h in hashes.items()])
 outputs=[p for p in o.iterdir() if p.suffix in {'.json','.tsv'}]
 dump(o/'validation.json',{'status':'PASS','checks':['input hashes unchanged and match original patient batch','full974 each model with674 planned','independent OLS residual rho for all computed rows','original effect/q unchanged','M0_CC exact primary reuse when inputs identical','model-specific complete-family BH','unavailable P/q remain missing'],'public_output_sha256':{p.name:sha(p) for p in outputs},'limitations':['Approximate Freedman-Lane exchangeability/homoskedasticity assumptions','No purity or technical batch adjustment','M3 unsupported sparse original categories; not mathematically rank deficient','Same CAMP patients, not independent validation']})
 (o/'.running').unlink();(o/'DONE').write_text('Completed scoped covariate follow-up\n');print(json.dumps(summary),flush=True)
if __name__=='__main__':
 try:main()
 except Exception:traceback.print_exc();raise
