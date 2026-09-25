"""SLC6A6 pooled-cell exploratory expression comparison; source cells stay server-side."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
from pathlib import Path
import sys,json,hashlib,importlib.util
import numpy as np,pandas as pd,h5py,scipy
from scipy import sparse,stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA')
def main():
 R=Path(sys.argv[1]);P=R/'public';P.mkdir();D=R/'private';D.mkdir();assert (R/'.running').exists()
 prev=BASE/'20260925T133935Z_lypla1_KEGGscore_v1';cf=prev/'source/brca_sc117_Wu2021_config.json';cfg=json.loads(cf.read_text());src=BASE/'20260919T140105Z_scRNA117_v1/source'/cfg['file']
 helper=prev/'brca_sc117_profile_v1.py';sp=importlib.util.spec_from_file_location('h',helper);h=importlib.util.module_from_spec(sp);sp.loader.exec_module(h)
 with h5py.File(src,'r') as f:
  obs=h.dataframe(f['obs']);var=h.dataframe(f['raw/var']);names=var.feature_name.astype(str);assert names.eq('SLC6A6').sum()==1;ix=int(np.flatnonzero(names.eq('SLC6A6'))[0]);groups=obs[cfg['celltype_column']].map(cfg['celltype_map']);sel=groups.isin(['Malignant_epithelial','Nonmalignant_epithelial']).to_numpy()
  raw=f['raw/X'];ptr=raw['indptr'][:];vals=[];labs=[]
  for start in range(0,len(obs),2000):
   end=min(start+2000,len(obs));mask=sel[start:end]
   if not mask.any():continue
   lo,hi=ptr[start],ptr[end];x=sparse.csr_matrix((raw['data'][lo:hi],raw['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(var)))[mask];lib=np.asarray(x.sum(1)).ravel().astype(float);assert (lib>0).all()
   vals.extend(np.log1p(x[:,ix].toarray().ravel().astype(float)/lib*10000));labs.extend(groups.iloc[start:end][mask])
 v=np.array(vals);m=np.array(labs)=='Malignant_epithelial';a,b=v[m],v[~m];assert len(a)==24489 and len(b)==4355
 np.savez_compressed(D/'Wu_cell_values.npz',expression=v,malignant=m)
 u,p=stats.mannwhitneyu(a,b,alternative='two-sided',method='asymptotic',use_continuity=True);eff=2*u/(len(a)*len(b))-1
 # Independent rank-sum construction and tie-corrected normal P.
 n,k=len(a),len(b);N=n+k;u2=stats.rankdata(v)[m].sum()-n*(n+1)/2;_,cnt=np.unique(v,return_counts=True);tie=np.sum(cnt.astype(float)**3-cnt);sd=np.sqrt(n*k/12*(N+1-tie/(N*(N-1))));z=(abs(u2-n*k/2)-.5)/sd;p2=2*stats.norm.sf(z);logp=(np.log(2)+stats.norm.logsf(z))/np.log(10)
 assert np.isclose(u,u2) and np.isclose(p,p2,rtol=1e-10,atol=1e-300)
 row=dict(cancer='BRCA',cohort='Wu2021',stage_id='06_EXTERNAL',run_id=R.name,analysis_version='slc6a6_cellcomparison_v1',analysis_type='pooled_cell_group_comparison',metabolite_key='NA',metabolite_name='NA',gene='SLC6A6',unit='cell',n=n,n_reference=k,effect_type='rank_biserial_malignant_minus_nonmalignant',effect=eff,ci_lower=np.nan,ci_upper=np.nan,p_value=p,q_value=p,test_family='SLC6A6_Wu_cell_level_single_test',family_n_evaluable=1,status='DONE',reason='EXPLORATORY_NO_DONOR_ADJUSTMENT',source_id=cfg['source_id'],mean_log1p_malignant=a.mean(),mean_log1p_nonmalignant=b.mean(),median_log1p_malignant=np.median(a),median_log1p_nonmalignant=np.median(b),mean_CP10k_malignant=np.expm1(a).mean(),mean_CP10k_nonmalignant=np.expm1(b).mean(),mean_CP10k_ratio=np.expm1(a).mean()/np.expm1(b).mean(),detection_malignant=(a>0).mean(),detection_nonmalignant=(b>0).mean(),log10_p=logp)
 save(pd.DataFrame([row]),P/'results.tsv')
 # Existing Pal source has no nonmalignant epithelial reference; do not fabricate one.
 cover=BASE/'20260925T140924Z_cellcomparison_v1/public/cohort_coverage.tsv';pd.read_csv(cover,sep='\t').to_csv(P/'cohort_coverage.tsv',sep='\t',index=False)
 save(pd.DataFrame([dict(cohort='Pal2021_reprocessed',gene='SLC6A6',status='NOT_EVALUABLE',reason='Current reprocessed source lacks annotated nonmalignant epithelial reference; no SLC6A6 comparison run')]),P/'not_evaluable.tsv')
 spec=dict(unit='pooled_cell',gene='SLC6A6',normalization='log1p(CP10k), float64',test='two-sided Mann-Whitney asymptotic, tie and continuity corrected',family_n=1,donor_adjustment=False,source_annotation=cfg['annotation_origin'],CI='NOT_RUN',software=dict(numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__),limitation='Cells within a donor are dependent; P is exploratory and can be anticonservative; no patient inference',parent_commit='85d7bbfc670911a47ee25e9a770e0d28355f20aa')
 (P/'analysis_spec.json').write_text(json.dumps(spec,indent=2));(P/'validation.json').write_text(json.dumps(dict(status='PASS',unique_exact_symbol=True,cell_counts_match_prior=True,independent_rank_U_and_P='PASS',figure_review='pending'),indent=2))
 manifest=[dict(path=str(q),sha256=sha(q)) for q in [src,cf,helper,cover,Path(__file__)]];save(pd.DataFrame(manifest),P/'source_manifest.tsv')
 plt.rcParams['pdf.fonttype']=42;fig,ax=plt.subplots(1,2,figsize=(9,4.5),constrained_layout=True);vp=ax[0].violinplot([b,a],showmedians=True,showextrema=False)
 for body,c in zip(vp['bodies'],['#4477AA','#EE6677']):body.set_facecolor(c);body.set_alpha(.75)
 ax[0].set_xticks([1,2],['Nonmalignant\nn='+str(k),'Malignant\nn='+str(n)]);ax[0].set_ylabel('SLC6A6 log1p(CP10k)');ax[0].set_title(f'Mean CP10k ratio = {row["mean_CP10k_ratio"]:.2f}\nCell-level P = {p:.3g}' if p else f'Mean CP10k ratio = {row["mean_CP10k_ratio"]:.2f}\nlog10(cell-level P) = {logp:.1f}')
 rates=[100*(b>0).mean(),100*(a>0).mean()];bars=ax[1].bar(['Nonmalignant','Malignant'],rates,color=['#4477AA','#EE6677']);ax[1].bar_label(bars,fmt='%.1f%%');ax[1].set_ylim(0,100);ax[1].set_ylabel('Cells with detected SLC6A6 (%)');ax[1].set_title('Detection fraction')
 fig.suptitle('Wu | SLC6A6 epithelial expression\nPooled cells; donor clustering not modelled',fontsize=12);fig.savefig(P/'SLC6A6_epithelial.png',dpi=200);fig.savefig(P/'SLC6A6_epithelial.pdf');plt.close(fig)
 print(json.dumps(row,indent=2),flush=True)
if __name__=='__main__':main()
