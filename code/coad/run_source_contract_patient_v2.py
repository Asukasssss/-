"""Server-only incremental patient analysis under the fixed cell-source contract."""
import argparse, json, hashlib, platform
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import scipy
from patient_statistics import association, bh, stable_seed

ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
VERSION='COAD_source_contract_patient_v2'
PREFIX='cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def rd(p):return pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False)
def save(t,p):
 t=t.copy(); cols=PREFIX+[x for x in t if x not in PREFIX]
 for c in PREFIX:
  if c not in t:t[c]='NA'
 t[cols].to_csv(p,sep='\t',index=False,na_rep='NA')
def key(r):return (r['metabolite_name'],r['metabolite_key'],r['human_gene_id'])
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();out=a.out
 assert out.parent==ROOT/'results/collaborative/COAD/B' and (out/'.running').exists()
 pub=out/'public';pub.mkdir(exist_ok=True);priv=out/'private';priv.mkdir(exist_ok=True)
 spec=json.loads((out/'locked_spec.json').read_text());assert spec['association_planned']==891 and spec['rna_planned']==526
 old=out.parent/'20260919T122034Z_patient_v1'; res=out.parent/'20260919T115842Z_identity_units_v1'
 manifest=rd(old/'source_manifest.tsv')
 for r in manifest.to_dict('records'):assert sha(ROOT/r['source_path'])==r['sha256'],r['source_path']
 clin=rd(res/'private_clinical_join.tsv');tm=clin[clin.TN.eq('Tumor')];nm=clin[clin.TN.eq('Normal')]
 people=tm[tm.stage.isin(['stage I','stage II','stage III','stage IV'])].copy(); normal=nm.set_index('individual').loc[people.individual]
 assert len(people)==33 and people.individual.is_unique and list(normal.index)==list(people.individual)
 base=ROOT/'data/candidates/camp_primary_tissue_multicancer'
 rp=base/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed/GSE89076.Agilent8x60K.log2_transformed.gene_symbol.csv'
 xp=base/'processed_metabolomics/PreprocessedData_COAD.xlsx'
 rna=pd.read_csv(rp,index_col=0); mats=pd.read_excel(xp,sheet_name=['metabo_imputed_filtered_Tumor','data'],index_col=0);met=mats['metabo_imputed_filtered_Tumor'];observed=mats['data']
 assert all(t.index.is_unique and t.columns.is_unique for t in [rna,met,observed])
 assert set(people.RNAID).issubset(rna.columns) and set(normal.RNAID).issubset(rna.columns)
 rel=rd(out/'main_relations.tsv'); membership=rd(out/'gene_membership.tsv'); assert len(rel)==891 and rel.human_gene_id.nunique()==526
 elig=rd(res/'candidate_eligibility.tsv');oldsymbols={r['human_gene_id']:r['rna_symbol'] for r in elig.to_dict('records')}
 def symbol(r):
  s=oldsymbols.get(r['human_gene_id'],'NA')
  if s!='NA':assert s in rna.index;return s
  return r['gene'] if r['gene'] in rna.index else 'NA'
 # Reuse is allowed only after ALL historical source hashes and public outputs match.
 oldcheck=json.loads((old/'validation.json').read_text())
 for name in ['results.tsv','paired_RNA33.tsv']:assert sha(old/name)==oldcheck['public_output_sha256'][name]
 prior={key(r):r for r in rd(old/'results.tsv').to_dict('records')}
 availpath=out.parent/'20260919T122940Z_availability_reuse_v2/results.tsv'
 assert sha(availpath)==spec['legacy_availability_sha256']
 priorav={key(r):r for r in rd(availpath).to_dict('records')}
 hashes={str(p):sha(p) for p in [rp,xp,res/'private_clinical_join.tsv',out/'main_relations.tsv',out/'gene_membership.tsv',out/'locked_spec.json',Path(__file__),out/'patient_statistics.py']}
 (pub/'analysis_spec.json').write_text(json.dumps(dict(spec,input_hashes=hashes,software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__)),indent=2))
 def base_row(r,kind):
  return dict(cancer='COAD',cohort='CAMP_COAD',stage_id='03_PATIENT',run_id=out.name,analysis_version=VERSION,analysis_type=kind,
   metabolite_key=r.get('metabolite_key','NA'),metabolite_name=r.get('metabolite_name','NA'),gene=r['gene'],human_gene_id=r['human_gene_id'],rna_symbol=symbol(r),
   unit='GEO_explicit_individual',effect_type='Spearman_rho',test_family=VERSION+'_'+kind,source_id=rp.name+';'+xp.name,
   relation_id='COAD|'+r.get('metabolite_name','NA')+'|'+r.get('metabolite_key','NA')+'|'+r['human_gene_id'],family_n_planned=891,reused_from='NA')
 rows=[]; selected_fields=['n','effect','p_value','ci_lower','ci_upper','bootstrap_valid','permutations','random_seed','n_unique_metabolite','n_unique_rna','loo_valid','loo_min','loo_max','loo_max_abs_delta','loo_any_sign_change']
 for i,r in enumerate(rel.to_dict('records')):
  s=symbol(r); main_result=None;main_good=None
  for kind,lookup in [('association_primary',prior),('association_available',priorav)]:
   row=base_row(r,kind)
   if s=='NA':row.update(status='NOT_EVALUABLE',reason='RNA_SYMBOL_UNAVAILABLE',n=0);rows.append(row);continue
   x=met.loc[r['metabolite_name'],people.MetabID].to_numpy(float);y=rna.loc[s,people.RNAID].to_numpy(float)
   if kind=='association_available':x=np.where(np.isfinite(observed.loc[r['metabolite_name'],people.MetabID].to_numpy(float)),x,np.nan)
   good=np.isfinite(x)&np.isfinite(y);vectorhash=hashlib.sha256(np.column_stack([x[good],y[good]]).astype('<f8').tobytes()+('|'.join(people.loc[good,'individual'])).encode()).hexdigest()
   cached=lookup.get(key(r))
   if cached and cached['status']=='DONE':
    assert cached['rna_symbol']==s and int(cached['n'])==int(good.sum())
    assert np.isclose(stats.spearmanr(x[good],y[good]).statistic,float(cached['effect']),atol=1e-12)
    row.update({f:cached[f] for f in selected_fields if f in cached});row.update(status='DONE',reason='REUSED_AFTER_SOURCE_HASH_AND_VECTOR_CHECK',reused_from='patient_v1' if kind=='association_primary' else 'availability_reuse_v2')
   elif kind=='association_available' and main_result is not None and np.array_equal(good,main_good):
    row.update({f:main_result[f] for f in selected_fields if f in main_result});row.update(status=main_result['status'],reason='IDENTICAL_INPUT_REUSED',reused_from='current_primary')
   else:
    row.update(association(x,y,stable_seed(VERSION+'|'+kind+'|'+row['relation_id']),loo=kind=='association_primary'))
   row['input_sha256']=vectorhash;rows.append(row)
   if kind=='association_primary':main_result=row.copy();main_good=good.copy()
  if (i+1)%100==0:print('RELATIONS',i+1,flush=True)
 ass=pd.DataFrame(rows)
 rr=[]
 for r in rel[['gene','human_gene_id']].drop_duplicates().to_dict('records'):
  row=base_row(r,'rna_paired_t');row.update(family_n_planned=526,effect_type='mean_paired_difference_log2_RNA',unit='explicit_individual_matched_pair')
  s=symbol(r)
  if s=='NA':row.update(status='NOT_EVALUABLE',reason='RNA_SYMBOL_UNAVAILABLE',n=0);rr.append(row);continue
  t=rna.loc[s,people.RNAID].to_numpy(float);n=rna.loc[s,normal.RNAID].to_numpy(float);good=np.isfinite(t)&np.isfinite(n);d=t[good]-n[good];nn=len(d)
  row.update(n=nn,n_reference=nn,n_up=int((d>0).sum()),n_down=int((d<0).sum()),n_equal=int((d==0).sum()),median_difference=float(np.median(d)) if nn else None)
  if nn<8 or np.std(d,ddof=1)==0:row.update(status='NOT_EVALUABLE',reason='INSUFFICIENT_PAIRS_OR_CONSTANT_DIFFERENCES');rr.append(row);continue
  mean=float(d.mean());se=float(stats.sem(d));lo,hi=stats.t.interval(.95,nn-1,loc=mean,scale=se);seed=stable_seed(VERSION+'|RNA|'+r['human_gene_id']);rng=np.random.default_rng(seed);boot=d[rng.integers(nn,size=(4000,nn))].mean(1)
  row.update(status='DONE',reason='NEW_METHOD_PAIRED_T_OLD_WILCOXON_PRESERVED',effect=mean,ci_lower=float(lo),ci_upper=float(hi),p_value=float(stats.ttest_1samp(d,0).pvalue),bootstrap_ci_lower=float(np.quantile(boot,.025)),bootstrap_ci_upper=float(np.quantile(boot,.975)),bootstrap_valid=4000,random_seed=seed,proportion_up=float((d>0).mean()),proportion_down=float((d<0).mean()))
  rr.append(row)
 rnas=pd.DataFrame(rr); summaries=[]
 for table in [ass,rnas]:
  table['q_value']=np.nan
  for kind,sub in table.groupby('analysis_type'):
   ix=sub.index[sub.status.eq('DONE')];pv=pd.to_numeric(table.loc[ix,'p_value']).to_numpy();q=bh(pv);table.loc[ix,'q_value']=q;table.loc[sub.index,'family_n_evaluable']=len(ix)
   from statsmodels.stats.multitest import multipletests
   assert np.allclose(q,multipletests(pv,method='fdr_bh')[1])
   summaries.append(dict(analysis_type=kind,planned=len(sub),evaluable=len(ix),p_lt_005=int((pv<.05).sum()),q_lt_005=int((q<.05).sum()),reused=int(sub.reused_from.ne('NA').sum())))
 save(ass,pub/'association_results.tsv');save(rnas,pub/'rna_results.tsv');save(rnas[pd.to_numeric(rnas.p_value,errors='coerce')<.05],pub/'rna_P_lt_005.tsv')
 save(ass[(pd.to_numeric(ass.p_value,errors='coerce')<.05)],pub/'association_P_lt_005.tsv')
 legacyrna=rd(old/'paired_RNA33.tsv');legacyrna.to_csv(pub/'historical_RNA_wilcoxon.tsv',sep='\t',index=False)
 coverage=membership.copy();coverage['rna_symbol']=[symbol(r) for r in coverage.to_dict('records')];coverage['rna_measurement_status']=np.where(coverage.rna_symbol.eq('NA'),'NOT_EVALUABLE','DONE');coverage.to_csv(pub/'gene_coverage.tsv',sep='\t',index=False)
 for p,h in hashes.items():assert sha(p)==h
 pd.DataFrame([dict(source_path=p,sha256=h) for p,h in hashes.items()]).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
 (pub/'summary.json').write_text(json.dumps(dict(status='DONE',families=summaries,patient_n=33,history_preserved=True),indent=2))
 (pub/'validation.json').write_text(json.dumps(dict(status='PASS',checks=['historical input hashes unchanged','historical output hashes unchanged','reuse rho and actual n checked against source vectors','BH independently checked','explicit patient pairing','unavailable values remain NA'],output_sha256={p.name:sha(p) for p in pub.iterdir() if p.is_file()}),indent=2))
 (out/'.running').unlink();(out/'DONE').write_text('DONE\n');print(json.dumps(summaries),flush=True)
if __name__=='__main__':main()
