"""Pooled-cell exploratory malignant versus nonmalignant epithelial comparison."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import sys,json,hashlib,importlib.util
from pathlib import Path
import numpy as np,pandas as pd,h5py,scipy
from scipy import sparse,stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A')
PREV=BASE/'20260925T133935Z_lypla1_KEGGscore_v1'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA')
def main():
 R=Path(sys.argv[1]);R.mkdir(exist_ok=False);(R/'.running').write_text('cellcomparison_v1')
 P=R/'public';P.mkdir();D=R/'private';D.mkdir()
 genes=json.loads((PREV/'public/membership_freeze.json').read_text())['genes'];assert len(genes)==92 and 'LYPLA1' not in genes
 spec=dict(version='cellcomparison_v1',unit='cell_not_donor',comparison='malignant vs annotated nonmalignant epithelial',gene_set='KEGG hsa00564 frozen92',genes=genes,normalization='per-cell log1p(CP10k)',score='mean gene-wise z across pooled malignant and nonmalignant cells within cohort; ddof1',test='two-sided Mann-Whitney asymptotic tie-corrected continuity-corrected; no donor adjustment',metrics=['LYPLA1','KEGG92_score'],selection='all annotated cells; no donor-count filter',ci='NOT_RUN: no donor-independent confidence intervals claimed',seed=20260925,software=dict(numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__),limitation='P values treat cells as independent and can be anticonservative; exploratory only, not patient inference')
 (P/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
 hs=importlib.util.spec_from_file_location('helper',PREV/'brca_sc117_profile_v1.py');h=importlib.util.module_from_spec(hs);hs.loader.exec_module(h)
 manifest=[dict(path=str(p),sha256=sha(p)) for p in [Path(__file__),PREV/'public/membership_freeze.json',PREV/'brca_sc117_profile_v1.py']]
 rows=[];coverage=[];fig,axes=plt.subplots(1,2,figsize=(11,4.6),constrained_layout=True);plt.rcParams['pdf.fonttype']=42
 for co,cn in [('Wu2021','brca_sc117_Wu2021_config.json'),('Pal2021_reprocessed','brca_sc117_Pal2021_config.json')]:
  cf=PREV/'source'/cn;cfg=json.loads(cf.read_text());src=BASE/'20260919T140105Z_scRNA117_v1/source'/cfg['file']
  manifest.extend([dict(path=str(cf),sha256=sha(cf)),dict(path=str(src),sha256=sha(src))])
  with h5py.File(src,'r') as f:
   obs=h.dataframe(f['obs']);var=h.dataframe(f[cfg['var_path']]);names=var[cfg['symbol_column']].astype(str)
   groups=obs[cfg['celltype_column']].map(cfg['celltype_map']);valid=np.ones(len(obs),bool)
   if cfg.get('filter_column'):valid &= obs[cfg['filter_column']].eq(cfg['filter_value']).to_numpy()
   counts={g:int((valid & groups.eq(g).to_numpy()).sum()) for g in ['Malignant_epithelial','Nonmalignant_epithelial']}
   coverage.append(dict(cohort=co,annotation_origin=cfg['annotation_origin'],**counts))
   base=dict(cancer='BRCA',cohort=co,stage_id='06_EXTERNAL',run_id=R.name,analysis_version='cellcomparison_v1',analysis_type='pooled_cell_group_comparison',metabolite_key='NA',metabolite_name='NA',unit='cell',n=counts['Malignant_epithelial'],n_reference=counts['Nonmalignant_epithelial'],effect_type='rank_biserial_malignant_minus_nonmalignant',ci_lower=np.nan,ci_upper=np.nan,test_family='two_metrics_Wu_pooled_cells',source_id=cfg['source_id'])
   if min(counts.values())==0:
    for metric in spec['metrics']:rows.append(dict(base,gene='LYPLA1' if metric=='LYPLA1' else 'NA',metric=metric,effect=np.nan,p_value=np.nan,q_value=np.nan,status='NOT_EVALUABLE',reason='No annotated nonmalignant epithelial reference in current Pal reprocessed source'))
    continue
   wanted=['LYPLA1']+genes;assert all(names.eq(g).sum()==1 for g in wanted)
   ix=[int(np.flatnonzero(names.eq(g))[0]) for g in wanted];sel=valid&groups.isin(counts).to_numpy();pieces=[];labels=[]
   raw=f[cfg['raw_path']];ptr=raw['indptr'][:]
   for start in range(0,len(obs),2000):
    end=min(start+2000,len(obs));mask=sel[start:end]
    if not mask.any():continue
    lo,hi=ptr[start],ptr[end];x=sparse.csr_matrix((raw['data'][lo:hi],raw['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(var)))[mask]
    lib=np.asarray(x.sum(1)).ravel();assert (lib>0).all()
    pieces.append(np.log1p(x[:,ix].toarray()/lib[:,None]*10000));labels.extend(groups.iloc[start:end][mask])
   X=np.vstack(pieces);lab=np.array(labels);m=lab=='Malignant_epithelial';sd=X[:,1:].std(0,ddof=1);assert (sd>0).all()
   score=((X[:,1:]-X[:,1:].mean(0))/sd).mean(1)
   np.savez_compressed(D/(co+'_cell_values.npz'),expression=X[:,0],score=score,malignant=m)
   save(pd.DataFrame(dict(gene=genes,mean=X[:,1:].mean(0),sd=sd)),P/'score_standardization.tsv')
   for j,(metric,v) in enumerate([('LYPLA1',X[:,0]),('KEGG92_score',score)]):
    a,b=v[m],v[~m];u,p=stats.mannwhitneyu(a,b,alternative='two-sided',method='asymptotic',use_continuity=True);eff=2*u/(len(a)*len(b))-1
    row=dict(base,gene='LYPLA1' if metric=='LYPLA1' else 'NA',metric=metric,effect=eff,p_value=p,q_value=np.nan,status='DONE',reason='CELL_LEVEL_EXPLORATORY_NO_DONOR_ADJUSTMENT',mean_malignant=a.mean(),mean_nonmalignant=b.mean(),median_malignant=np.median(a),median_nonmalignant=np.median(b),mean_difference=a.mean()-b.mean(),detection_malignant=(a>0).mean() if metric=='LYPLA1' else np.nan,detection_nonmalignant=(b>0).mean() if metric=='LYPLA1' else np.nan,mean_CP10k_ratio=np.expm1(a).mean()/np.expm1(b).mean() if metric=='LYPLA1' else np.nan,p_underflow=bool(p==0))
    rows.append(row);ax=axes[j];parts=ax.violinplot([b,a],showmedians=True,showextrema=False,points=100)
    for body,c in zip(parts['bodies'],['#4477AA','#EE6677']):body.set_facecolor(c);body.set_alpha(.75)
    ax.set_xticks([1,2],['Nonmalignant\nn='+str(len(b)),'Malignant\nn='+str(len(a))]);ax.set_ylabel('log1p(CP10k)' if j==0 else 'Mean per-cell gene z-score');ps='below floating-point precision' if p==0 else f'{p:.3g}'
    ax.set_title(f'Wu | {metric}\nrank-biserial={eff:.3f}; cell-level P={ps}',fontsize=10)
   print(co,counts,flush=True)
 df=pd.DataFrame(rows);ok=df.status.eq('DONE');pv=df.loc[ok,'p_value'].to_numpy();order=np.argsort(pv);q=np.empty(len(pv));q[order]=np.minimum.accumulate((pv[order]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1].clip(0,1);df.loc[ok,'q_value']=q;df['family_n_evaluable']=len(pv)
 save(df,P/'results.tsv');save(pd.DataFrame(coverage),P/'cohort_coverage.tsv');save(pd.DataFrame(manifest),P/'source_manifest.tsv')
 fig.suptitle('Pooled epithelial cells: exploratory; donor clustering NOT modelled',fontsize=12)
 fig.savefig(P/'malignant_vs_nonmalignant.png',dpi=200);fig.savefig(P/'malignant_vs_nonmalignant.pdf');plt.close(fig)
 assert len(df)==4 and len(pv)==2 and np.isfinite(df.loc[ok,'effect']).all()
 (P/'validation.json').write_text(json.dumps(dict(status='PASS',n_tests=2,n_not_evaluable=2,unique_exact_gene_symbols=True,frozen92=True,no_donor_inference=True,independent_recalculation='pending'),indent=2))
 print(df[['cohort','metric','n','n_reference','effect','p_value','mean_CP10k_ratio']].to_string(index=False),flush=True)
if __name__=='__main__':main()
