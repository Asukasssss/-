"""Focused LYPLA1 comparison using author patient labels, never cells as replicates."""
import argparse, hashlib, json, itertools, platform
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
from scipy import sparse,stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''): h.update(b)
 return h.hexdigest()
def main(out,commit):
 out.mkdir(parents=True,exist_ok=True)
 with (out/'.running').open('x') as f:f.write('LYPLA1')
 pub=out/'public';pub.mkdir(exist_ok=False);priv=out/'private';priv.mkdir(exist_ok=False)
 s=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/gse159115_kyn_cell_source_v0.1')
 aps=[s/'GSE159115_ccRCC_anno.csv.gz',s/'GSE159115_normal_anno.csv.gz']
 a=pd.read_csv(aps[0]);b=pd.read_csv(aps[1]);a=a[a.anno.eq('Tumor')].copy();a['tissue']='Tumor';b['tissue']='Normal'
 obs=pd.concat([a,b]).set_index('cell');assert obs.index.is_unique and obs.patient.notna().all()
 obs['expr']=np.nan;obs['det']=np.nan;used=[];seen=set()
 pt=['PT-A','PT-B','PT-C'];tub=pt+['TAL','tAL','DCT','CNT','DL','PC','IC-A','IC-B','IC-PC']
 spec=dict(gene='LYPLA1',stable_gene_id='ENSG00000120992',code_commit=commit,version='ccrcc_lypla1_epithelial_v1',primary='Tumor vs pooled author PT-A/PT-B/PT-C; matched author patient labels',sensitivity='Tumor vs pooled defined renal tubular epithelium',normal_tubular_labels=tub,exclusions='UC,Podo,unknown,ua and non-tubular categories excluded from broad comparator',minimum_cells_per_donor_group=20,normalization='log1p(10000*counts/full_H5_library)',unit='patient',test='two-sided exhaustive paired sign permutation on mean patient difference',family='two contrasts; BH over two P values',ci='95% paired t interval; small-n normality assumption',all_donor_summaries='descriptive only, overlapping donors not treated as independent groups',seed=20260926,source_url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159115',software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__))
 (pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
 for p in sorted((s/'h5').glob('*.h5')):
  sample='SI_'+p.name.split('_SI_',1)[1].split('_filtered')[0]
  if sample not in set(obs['sample']):continue
  with h5py.File(p,'r') as f:
   z=f['GRCh38'];ids=[v.decode().split('.')[0] for v in z['genes'][:]];ii=[i for i,x in enumerate(ids) if x==spec['stable_gene_id']];assert len(ii)==1
   bars=[sample+'_'+v.decode() for v in z['barcodes'][:]];x=sparse.csc_matrix((z['data'][:],z['indices'][:],z['indptr'][:]),shape=(len(ids),len(bars)))
   j=np.flatnonzero(pd.Index(bars).isin(obs.index));cells=np.array(bars)[j];assert not seen.intersection(cells);seen.update(cells);x=x[:,j]
   assert np.isfinite(x.data).all() and (x.data>=0).all() and (x.data==np.floor(x.data)).all()
   lib=np.asarray(x.sum(axis=0)).ravel();assert (lib>0).all();v=x[ii[0]].toarray().ravel()
   obs.loc[cells,'expr']=np.log1p(v*10000/lib);obs.loc[cells,'det']=(v>0).astype(float)
  used.append(p);print('read',p.name,flush=True)
 assert seen==set(obs.index) and obs.expr.notna().all()
 profiles=[]
 groups={'Tumor':obs.tissue.eq('Tumor'),'Normal_PT':obs.tissue.eq('Normal')&obs.anno.isin(pt),'Normal_tubular':obs.tissue.eq('Normal')&obs.anno.isin(tub)}
 for c in pt:groups['Normal_'+c]=obs.tissue.eq('Normal')&obs.anno.eq(c)
 for label,mask in groups.items():
  for donor,z in obs[mask].groupby('patient'):
   profiles.append(dict(patient=donor,group=label,cells=len(z),expression=z.expr.mean(),detection=z.det.mean()))
 pro=pd.DataFrame(profiles);pro.to_csv(priv/'donor_profiles.tsv',sep='\t',index=False)
 eligible=pro[pro.cells.ge(20)];summ=[];rows=[];fig,axes=plt.subplots(1,2,figsize=(9,4))
 for label,z in eligible.groupby('group'):
  summ.append(dict(group=label,n_donors=len(z),n_cells=z.cells.sum(),mean_expression=z.expression.mean(),mean_detection=z.detection.mean(),status='DONE' if len(z)>=3 else 'DESCRIPTIVE_LOW_N'))
 for ax,normal in zip(axes,['Normal_PT','Normal_tubular']):
  t=eligible[eligible.group.eq('Tumor')].set_index('patient');n=eligible[eligible.group.eq(normal)].set_index('patient');shared=sorted(set(t.index)&set(n.index));tv=t.loc[shared].expression.to_numpy();nv=n.loc[shared].expression.to_numpy();d=tv-nv;k=len(d)
  perm=np.array(list(itertools.product([-1,1],repeat=k)));p=float(np.mean(np.abs((perm*d).mean(1))>=abs(d.mean())-1e-14));se=stats.sem(d);ci=stats.t.interval(.95,k-1,loc=d.mean(),scale=se)
  rows.append(dict(cancer='ccRCC',cohort='GSE159115',stage_id='06_EXTERNAL',run_id=out.name,analysis_version=spec['version'],analysis_type='paired_epithelial_expression',metabolite_key='NA',metabolite_name='NA',gene='LYPLA1',unit='author_patient',n=k,n_reference=k,effect_type='Tumor-minus-normal mean log1p_CP10K',effect=d.mean(),ci_lower=ci[0],ci_upper=ci[1],p_value=p,q_value=np.nan,test_family='LYPLA1_two_epithelial_contrasts',family_n_evaluable=2,status='DONE',reason='small_n_exploratory',source_id='GSE159115',comparator=normal,tumor_mean=tv.mean(),normal_mean=nv.mean(),tumor_detection=t.loc[shared].detection.mean(),normal_detection=n.loc[shared].detection.mean(),tumor_cells=t.loc[shared].cells.sum(),normal_cells=n.loc[shared].cells.sum(),n_up=int((d>0).sum()),n_down=int((d<0).sum())))
  for x,y in zip(nv,tv):ax.plot([0,1],[x,y],'-o',color='#627d98',alpha=.7)
  ax.set_xticks([0,1]);ax.set_xticklabels([normal,'Tumor']);ax.set_ylabel('Patient mean log1p(CP10K)');ax.set_title('LYPLA1 | paired n=%d | P=%.3g'%(k,p));ax.spines[['top','right']].set_visible(False)
 df=pd.DataFrame(rows);p=df.p_value.to_numpy();order=np.argsort(p);q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];df.loc[order,'q_value']=np.minimum(q,1)
 df.to_csv(pub/'results.tsv',sep='\t',index=False);pd.DataFrame(summ).to_csv(pub/'group_summaries.tsv',sep='\t',index=False)
 # Plot is aggregate visualization; no donor identifiers or per-donor table exported.
 fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(pub/('LYPLA1_paired_epithelium.'+ext),dpi=180)
 plt.close(fig)
 manifest=[dict(source_id=p.name,path_or_url=str(p),sha256=sha(p),access_scope='SERVER_PRIVATE' if p in aps+used else 'CODE') for p in aps+used+[Path(__file__)]]
 pd.DataFrame(manifest).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
 (pub/'validation.json').write_text(json.dumps(dict(status='DONE',all_cells_matched=True,integer_counts=True,matched_by='author patient field',raw_cell_overlap=0,contrast_count=len(df),p_values_exact=True,donor_table_server_only=True),indent=2))
 (out/'.running').rename(out/'.done');print(df.to_string(index=False),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--commit',required=True);a=p.parse_args();main(a.out,a.commit)
