"""Server-only incremental CAMP analysis; old statistics remain immutable."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
from pathlib import Path
import sys,json,hashlib,platform
import numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from camp_per_cancer_repro import calculate
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
BASE=ROOT/'results/collaborative/BRCA/A'
R=Path(sys.argv[1]); P=R/'public'; V='mapping262_patient_v1'
PREFIX='cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()
def save(d,p): d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def seed(s): return int(hashlib.sha256((V+'|'+s).encode()).hexdigest()[:8],16)
def base(g,f,k='NA',name='NA'):
 d=dict.fromkeys(PREFIX,np.nan);d.update(cancer='BRCA',cohort='CAMP_BRCA1_Terunuma',stage_id='03_PATIENT',run_id=R.name,analysis_version=V,analysis_type=f,metabolite_key=k,metabolite_name=name,gene=g,unit='author_case_row',effect_type='Spearman_rho',test_family=f,status='NOT_EVALUABLE',reason='feature_missing',source_id='CAMP_GSE37751_original_author_case_table',relation_id=k+'|'+g,statistics_reused=False)
 return d
def correct(d):
 for f,ix in d.groupby('test_family').groups.items():
  ok=d.loc[ix,'p_value'].notna(); ii=np.array(list(ix))[ok];d.loc[ix,'family_n_evaluable']=len(ii)
  if len(ii):
   ps=d.loc[ii,'p_value'].to_numpy();q=multipletests(ps,method='fdr_bh')[1];order=np.argsort(ps);check=np.minimum(1,np.minimum.accumulate((ps[order]*len(ps)/np.arange(1,len(ps)+1))[::-1])[::-1]);assert np.allclose(q[order],check,rtol=0,atol=1e-12);d.loc[ii,'q_value']=q
 return d[PREFIX+[c for c in d if c not in PREFIX]]
def main():
 assert (R/'.running').read_text().strip()==V and not (P/'CAMP_relations262.tsv').exists()
 rel=pd.read_csv(R/'source/direct_relations.tsv',sep='\t');assert len(rel)==262 and rel.relation_key.is_unique and rel.gene.nunique()==150
 olddir=BASE/'20260921T103911Z_camp_pair_sensitivity_v1/public'
 old=pd.read_csv(olddir/'association_174_sensitivities.tsv',sep='\t').set_index(['relation_id','test_family']);oe=pd.read_csv(olddir/'paired_RNA117.tsv',sep='\t').set_index('gene')
 identity=BASE/'20260921T102429Z_camp_sample_identity_v1/private/audited_mapping.tsv';m=pd.read_csv(identity,sep='\t',dtype=str);bad=m.TN.ne(m.pdf_TN);assert bad.sum()==1 and m.loc[bad,'GSM'].iloc[0]=='GSM927051'
 tm=m[m.TN.eq('Tumor')&~bad];nm=m[m.TN.eq('Normal')&~bad];assert len(tm)==60 and tm.case_row.is_unique and nm.case_row.is_unique
 pairs=tm[['case_row','RNAID']].merge(nm[['case_row','RNAID']],on='case_row',suffixes=('_t','_n'),validate='one_to_one');assert len(pairs)==45
 src=ROOT/'data/candidates/camp_primary_tissue_multicancer';rf=src/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/m.RNAFile.iloc[0];mf=src/'processed_metabolomics'/m.MetabFile.iloc[0]
 rna=pd.read_csv(rf,index_col=0);met=pd.read_excel(mf,sheet_name=tm.MetabFile_sheet.iloc[0],index_col=0);obs=pd.read_excel(mf,sheet_name='data',index_col=0)
 for d in [rna,met,obs]:d.columns=d.columns.astype(str);assert d.index.is_unique and d.columns.is_unique
 for d in [met,obs]:d.index=d.index.astype(str).str.strip();assert d.index.is_unique
 # Historical LPCAT4 can denote MBOAT2; do not infer modern Q643R3 identity from a symbol.
 blocked={'LPCAT4':'historical_LPCAT4_symbol_ambiguous_Q643R3_vs_Q6ZWT7;stable_probe_identity_required'}
 spec=dict(version=V,source_commit='9e4387f',relations=262,genes=150,history_union=156,primary='60 author tumor cases; documented tissue conflict excluded',paired_RNA='45 explicit author-case tumor-normal pairs; unchanged expression scale; paired t test',families=['processed_Spearman262','available_Spearman262','paired_RNA150'],BH='all evaluable of262 separately per correlation family; all evaluable of150 for RNA',reuse='same original matrix files, sample mapping, exact relation and method; n/effect checked; original effect/P/CI copied; new q',new_Spearman='same helper:9999 permutations plus1;4000 case-bootstrap percentile95CI;min8',paired_bootstrap=4000,metabolite_name_normalization='trim outer whitespace only; unique index checked',gene_alias_policy='no unsupported alias substitution; LPCAT4 blocked pending stable probe identity',postselection_exploratory=True,patient_values_exported=False,python=platform.python_version())
 (P/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
 paths=[identity,rf,mf,olddir/'association_174_sensitivities.tsv',olddir/'paired_RNA117.tsv',R/'source/direct_relations.tsv',Path(__file__),R/'camp_per_cancer_repro.py'];hashes={str(p):sha(p) for p in paths}
 cov=[];rows=[]
 for i,t in enumerate(rel.itertuples()):
  name=t.metabolite_name.strip();good=t.gene in rna.index and name in met.index and name in obs.index and t.gene not in blocked
  if good:x=met.loc[name,tm.MetabID].to_numpy(float);y=rna.loc[t.gene,tm.RNAID].to_numpy(float);avail=np.isfinite(obs.loc[name,tm.MetabID].to_numpy(float))
  for fam in ['processed_Spearman','available_Spearman']:
   b=base(t.gene,fam+'262',t.metabolite_key,t.metabolite_name);b['in_original117']=t.gene_in_original117
   pr=old.loc[(t.relation_key,fam)] if (t.relation_key,fam) in old.index else None
   if pr is not None:b.update(previous_p=pr.p_value,previous_q=pr.q_value,previous_effect=pr.effect,previous_run=olddir.parent.name)
   if t.gene in blocked:b.update(status='NEEDS_REVIEW',reason=blocked[t.gene])
   elif good:
    mask=np.isfinite(x)&np.isfinite(y)&(avail if fam.startswith('available') else True);xx=x[mask];yy=y[mask];b['n']=len(xx)
    if len(xx)>=8 and np.ptp(xx)>0 and np.ptp(yy)>0:
     rho=stats.spearmanr(xx,yy).statistic
     if pr is not None and pd.notna(pr.p_value):
      assert len(xx)==pr.n and abs(rho-pr.effect)<1e-12,(t.relation_key,fam)
      b.update({c:pr[c] for c in ['effect','ci_lower','ci_upper','p_value']});b.update(statistics_reused=True,seed=pr.get('seed',np.nan))
     else:
      sd=seed(fam+'|'+t.relation_key);z=calculate(xx,yy,sd);b.update(effect=z['rho'],p_value=z['p_value'],ci_lower=z['ci_lower'],ci_upper=z['ci_upper'],seed=sd,bootstrap_valid=z['bootstrap_valid'])
     b.update(status='DONE',reason='same_cohort_postselection;RNA_not_flux')
    else:b['reason']='less_than8_or_constant'
   rows.append(b)
  if i%40==0:print('CAMP',i+1,flush=True)
 save(correct(pd.DataFrame(rows)),P/'CAMP_relations262.tsv')
 er=[]
 for g in sorted(rel.gene.unique()):
  b=base(g,'paired_RNA150');b.update(effect_type='mean_tumor_minus_normal_author_expression_scale',reason='RNA_gene_missing')
  pr=oe.loc[g] if g in oe.index else None
  if pr is not None:b.update(previous_p=pr.p_value,previous_q=pr.q_value,previous_effect=pr.effect,previous_run=olddir.parent.name)
  if g in blocked:b.update(status='NEEDS_REVIEW',reason=blocked[g])
  elif g in rna.index:
   x=rna.loc[g,pairs.RNAID_t].to_numpy(float);y=rna.loc[g,pairs.RNAID_n].to_numpy(float);mask=np.isfinite(x)&np.isfinite(y);delta=x[mask]-y[mask];n=len(delta);b.update(n=n,n_reference=n)
   if n>=8 and delta.std(ddof=1)>0:
    mean=float(delta.mean());test=stats.ttest_rel(x[mask],y[mask]);se=delta.std(ddof=1)/np.sqrt(n);ci=stats.t.interval(.95,n-1,loc=mean,scale=se)
    if pr is not None and pd.notna(pr.p_value):
     assert pr.n==n and abs(pr.effect-mean)<1e-12 and np.isclose(pr.p_value,test.pvalue,rtol=1e-10,atol=1e-15)
     b.update({c:pr[c] for c in ['effect','ci_lower','ci_upper','p_value','bootstrap_mean_lower','bootstrap_mean_upper']});b['statistics_reused']=True
    else:
     sd=seed('paired_RNA|'+g);rng=np.random.default_rng(sd);bs=delta[rng.integers(n,size=(4000,n))].mean(1);b.update(effect=mean,p_value=test.pvalue,ci_lower=ci[0],ci_upper=ci[1],bootstrap_mean_lower=np.quantile(bs,.025),bootstrap_mean_upper=np.quantile(bs,.975),seed=sd)
    b.update(status='DONE',reason='paired_t_on_author_scale;not_protein_or_flux',positive_pairs=int((delta>0).sum()),negative_pairs=int((delta<0).sum()),equal_pairs=int((delta==0).sum()),positive_pair_fraction=float((delta>0).mean()),median_paired_difference=float(np.median(delta)))
   else:b['reason']='less_than8_or_constant'
  er.append(b)
 e=correct(pd.DataFrame(er));save(e,P/'paired_RNA150.tsv')
 save(pd.DataFrame([dict(gene=g,exact_symbol_in_CAMP=g in rna.index,status='NEEDS_REVIEW' if g in blocked else ('DONE' if g in rna.index else 'NOT_EVALUABLE'),reason=blocked.get(g,'exact_symbol_only_no_alias_expansion')) for g in sorted(rel.gene.unique())]),P/'gene_identity_coverage.tsv')
 assert all(sha(p)==h for p,h in hashes.items());save(pd.DataFrame([dict(path=p,sha256=h) for p,h in hashes.items()]),P/'source_manifest.tsv')
 a=pd.concat([correct(pd.DataFrame(rows)),e],ignore_index=True);v=dict(status='DONE',source_hashes_unchanged=True,old_statistics_unchanged=True,old_evaluable_overlap_checked=True,BH_independently_checked=True,tumor_cases=60,paired_cases=45,families={f:dict(planned=len(z),evaluable=int(z.p_value.notna().sum()),reused=int(z.statistics_reused.sum()),newly_calculated=int((z.p_value.notna()&~z.statistics_reused).sum()),q_lt005=int(z.q_value.lt(.05).sum())) for f,z in a.groupby('test_family')})
 (P/'validation.json').write_text(json.dumps(v,indent=2));(R/'NUMERICAL_DONE').write_text('DONE');print(json.dumps(v),flush=True)
if __name__=='__main__':main()
