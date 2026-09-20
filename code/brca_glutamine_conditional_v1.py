"""Execute frozen six-model exploratory partial-rank analysis on server-only input."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import argparse,json,hashlib,sys,itertools
from pathlib import Path
import numpy as np,pandas as pd,scipy
from scipy.stats import rankdata,t as tdist,spearmanr
from statsmodels.stats.multitest import multipletests

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(df,p):df.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def partial(a):
 r=np.column_stack([rankdata(a[:,i]) for i in range(a.shape[1])]);z=np.column_stack([np.ones(len(r)),r[:,2:]])
 if np.linalg.matrix_rank(z)!=z.shape[1]:raise ValueError('rank_deficient_design')
 e=r[:,:2]-z@np.linalg.lstsq(z,r[:,:2],rcond=None)[0]
 if np.min(np.std(e,axis=0))<1e-10:raise ValueError('constant_residual')
 rho=float(np.corrcoef(e.T)[0,1])
 # Independent precision-matrix identity for partial correlation.
 inv=np.linalg.inv(np.corrcoef(r.T));check=-inv[0,1]/np.sqrt(inv[0,0]*inv[1,1])
 assert abs(rho-check)<1e-10
 return rho

def main(root):
 pub=root/'public';pub.mkdir(exist_ok=True)
 spec=json.loads((root/'analysis_spec.json').read_text());src=Path(spec['server_input']);d=pd.read_csv(src,sep='\t')
 assert d.patient_id.is_unique
 names=spec['genes']+[spec['external_peak']];assert set(names)<=set(d)
 a=d[['patient_id']+names].replace([np.inf,-np.inf],np.nan).dropna();assert len(a)>=20
 out=[];checks=[];rng=np.random.default_rng(spec['seed'])
 for m in spec['models']:
  fields=[m['gene'],spec['external_peak']]+m['covariates'];x=a[fields].to_numpy(float);k=len(m['covariates']);rho=partial(x)
  df=len(x)-k-2;p=2*tdist.sf(abs(rho)*np.sqrt(df/(1-rho*rho)),df)
  boot=[]
  for i in range(2000):
   try:boot.append(partial(x[rng.integers(len(x),size=len(x))]))
   except (ValueError,np.linalg.LinAlgError):pass
  assert len(boot)>=1900
  lo,hi=np.quantile(boot,[.025,.975])
  base=float(spearmanr(x[:,0],x[:,1]).statistic)
  r=dict(cancer='BRCA',cohort='FUSCC_TNBC',stage_id='04_ROBUSTNESS',run_id=root.name,analysis_version='glutamine_context_v1',analysis_type=m['model'],metabolite_key=spec['metabolite_key'],metabolite_name='glutamine',gene=m['gene'],unit='author_patient_id',n=len(x),n_reference=np.nan,effect_type='partial_rank_rho',effect=rho,ci_lower=lo,ci_upper=hi,p_value=p,q_value=np.nan,test_family='exploratory_conditional6',family_n_evaluable=6,status='DONE',reason='post_selection_exploratory;approximate_t_p;pointwise_bootstrap_CI',source_id='FUSCC_BRCA_2022',covariates=';'.join(m['covariates']),baseline_same_subset_rho=base,df=df,bootstrap_valid=len(boot))
  out.append(r)
  b=a[~a.patient_id.isin(spec['warning_ids'])][fields].to_numpy(float)
  checks.append(dict(gene=m['gene'],model=m['model'],n=len(b),primary_rho=rho,exclude_warning_rho=partial(b),p_value='NA',q_value='NA',reason='prespecified_effect_only_sensitivity'))
  print(m['gene'],m['model'],rho,flush=True)
 tab=pd.DataFrame(out);tab['q_value']=multipletests(tab.p_value,method='fdr_bh')[1];tab['holm_p']=multipletests(tab.p_value,method='holm')[1]
 save(tab,pub/'conditional6_results.tsv');save(pd.DataFrame(checks),pub/'warning_sensitivity.tsv')
 pairs=[dict(gene1=g,gene2=h,n=len(a),rho=spearmanr(a[g],a[h]).statistic,interpretation='descriptive_expression_correlation;no_new_hypothesis_test') for g,h in itertools.combinations(spec['genes'],2)]
 save(pd.DataFrame(pairs),pub/'gene_pair_descriptive.tsv')
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,ax=plt.subplots(figsize=(8,4.5));y=np.arange(len(tab));ax.errorbar(tab.effect,y,xerr=np.array([tab.effect-tab.ci_lower,tab.ci_upper-tab.effect]),fmt='o',color='#246b91',capsize=3)
 ax.set_yticks(y);ax.set_yticklabels([f"{r.gene} | {r.covariates}" for r in tab.itertuples()]);ax.invert_yaxis();ax.axvline(0,color='gray',lw=1);ax.set_xlabel('Partial rank correlation (pointwise bootstrap 95% CI)');ax.set_title('FUSCC glutamine | exploratory conditional analysis');ax.spines[['top','right']].set_visible(False);fig.tight_layout();fig.savefig(pub/'conditional6_forest.png',dpi=180);plt.close(fig)
 val=dict(status='DONE',n_input=len(d),n_complete=len(a),n_sensitivity=int((~a.patient_id.isin(spec['warning_ids'])).sum()),models=6,bootstrap_resamples_per_model=2000,precision_matrix_crosscheck='PASS',q_family=6,prior_statistics_modified=False,patient_values_exported=False,independent_replication=False,python=sys.version,numpy=np.__version__,scipy=scipy.__version__)
 (pub/'validation.json').write_text(json.dumps(val,indent=2)+'\n')
 save(pd.DataFrame([dict(path=str(p),sha256=sha(p)) for p in [src,root/'analysis_spec.json',Path(__file__)]]),pub/'source_manifest.tsv')
 print(json.dumps(val),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args().root)
