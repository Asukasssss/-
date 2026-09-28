"""Author-labeled malignant epithelium LYPLA1 analysis; run only on server165."""
from pathlib import Path
import json,hashlib,platform
import h5py,numpy as np,pandas as pd,scipy
from scipy import sparse,stats
RUN=Path(__file__).resolve().parent;PUB=RUN/'public';PRIV=RUN/'private'
INPUT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/SCP1089/GSE202051_totaldata-final-toshare.h5ad')
SEED=20260928
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def decode(v):return v.decode() if isinstance(v,bytes) else str(v)
def read(o,k):
 v=o[k][:]
 if '__categories' in o and k in o['__categories']:
  cat=[decode(x) for x in o['__categories'][k][:]];return np.array([cat[i] if i>=0 else None for i in v],object)
 if v.dtype.kind in ['O','S','U']:return np.array([decode(x) for x in v],object)
 return v
def bh(p):
 p=np.asarray(p,float);q=np.full(len(p),np.nan);ix=np.where(np.isfinite(p))[0];ix=ix[np.argsort(p[ix])]
 if len(ix):q[ix]=np.minimum.accumulate((p[ix]*len(ix)/np.arange(1,len(ix)+1))[::-1])[::-1].clip(0,1)
 return q
def write(df,name):df.to_csv(PUB/name,sep='\t',index=False,na_rep='NA')
def extract():
 with h5py.File(INPUT,'r') as f:
  keys=['_index','pid','sampleid','new_treatment','Level 1 Annotation','Level 2 Annotation','Level 3 Annotation','total_counts','n_counts','cnv_score']
  c=pd.DataFrame({k:read(f['obs'],k) for k in keys}).rename(columns={'_index':'cell','new_treatment':'treatment','Level 1 Annotation':'level1','Level 2 Annotation':'level2','Level 3 Annotation':'level3'})
  assert len(c)==224988 and c.cell.is_unique and c.pid.notna().all() and (c.groupby('pid').treatment.nunique()==1).all()
  genes=read(f['var'],'_index');gi=np.where(genes=='LYPLA1')[0];assert len(gi)==1;gi=gi[0]
  cnt=f['layers/counts'];xx=f['X'];ptr=cnt['indptr'][:];assert np.array_equal(ptr,xx['indptr'][:])
  assert tuple(cnt.attrs['shape'])==(len(c),len(genes))
  values=np.zeros(len(c));expression=np.zeros(len(c));library=np.zeros(len(c));seen=np.zeros(len(c),bool)
  for lo in range(0,len(c),5000):
   hi=min(lo+5000,len(c));start,end=int(ptr[lo]),int(ptr[hi]);idx=cnt['indices'][start:end];v=cnt['data'][start:end]
   assert np.isfinite(v).all() and (v>=0).all() and np.equal(v,np.floor(v)).all()
   assert np.array_equal(idx,xx['indices'][start:end]);x=xx['data'][start:end];assert np.isfinite(x).all() and (x>=0).all()
   mat=sparse.csr_matrix((v,idx,ptr[lo:hi+1]-start),shape=(hi-lo,len(genes)))
   library[lo:hi]=np.asarray(mat.sum(axis=1)).ravel()
   hits=np.where(idx==gi)[0];rows=np.searchsorted(ptr,hits+start,side='right')-1
   assert len(np.unique(rows))==len(rows) and not seen[rows].any();seen[rows]=True
   values[rows]=v[hits];expression[rows]=x[hits]
   if lo%25000==0:print('EXTRACT',hi,flush=True)
  assert (library>0).all();c['gene_count']=values;c['matrix_library']=library;c['expression']=expression;c['detected']=values>0
  checks={}
  for name,lib in [('matrix_library',library),('total_counts',c.total_counts.to_numpy()),('n_counts',c.n_counts.to_numpy())]:
   good=lib>0;pred=np.log1p(10000*values[good]/lib[good]);checks[name]=bool(np.allclose(pred,expression[good],rtol=3e-5,atol=3e-6))
  assert any(checks.values()),'Author expression normalization not verified'
  assert ((expression>0)==(values>0)).all()
  assert (c.loc[c.level2=='Malignant','level1']=='Epithelial (malignant)').all()
  assert (c.loc[c.level2.isin(['Ductal','Ductal (atypical)']),'level1']=='Epithelial (non-malignant)').all()
 c.to_csv(PRIV/'cell_values.tsv.gz',sep='\t',index=False)
 (PRIV/'normalization_checks.json').write_text(json.dumps(checks));print('NORMALIZATION',checks,flush=True)
 return c,checks
def paired(x,y,seed):
 d=x-y;rng=np.random.default_rng(seed);extreme=0
 for _ in range(100):
  perm=(rng.choice([-1.,1.],size=(1000,len(d)))*d).mean(axis=1);extreme+=int((abs(perm)>=abs(d.mean())-1e-12).sum())
 boot=d[rng.integers(0,len(d),(20000,len(d)))].mean(axis=1);lo,hi=np.quantile(boot,[.025,.975])
 return dict(effect=d.mean(),p_value=(extreme+1)/100001,ci_lower=lo,ci_upper=hi,higher_donors=int((d>0).sum()),lower_donors=int((d<0).sum()))
def analyze(c):
 rows=[];profiles=[];hist=[];donor=[]
 epi=c[c.level2.isin(['Malignant','Ductal','Ductal (atypical)'])]
 for (tr,ct),b in epi.groupby(['treatment','level2']):
  for pop in ['All_nuclei','LYPLA1_positive_nuclei']:
   d=b if pop=='All_nuclei' else b[b.detected];x=d.expression.to_numpy();u=d.groupby('pid').expression.agg(['mean','size']);u=u[u['size']>=20]
   rng=np.random.default_rng(SEED);uv=u['mean'].to_numpy();ci=np.quantile(uv[rng.integers(0,len(uv),(20000,len(uv)))].mean(axis=1),[.025,.975]) if len(uv)>=5 else [np.nan,np.nan]
   profiles.append(dict(treatment=tr,celltype=ct,population=pop,n_cells=len(d),n_donors=d.pid.nunique(),n_samples=d.sampleid.nunique(),
     n_eligible_donors=len(u),pooled_mean=x.mean() if len(x) else np.nan,pooled_median=np.median(x) if len(x) else np.nan,
     detection_fraction=b.detected.mean(),all_nuclei=len(b),positive_nuclei=int(b.detected.sum()),
     donor_mean=uv.mean() if len(uv) else np.nan,donor_ci_lower=ci[0],donor_ci_upper=ci[1],largest_donor_fraction=d.pid.value_counts().max()/len(d) if len(d) else np.nan))
   edges=np.linspace(0,10,401);counts,_=np.histogram(x,bins=edges);assert counts.sum()==len(d),'Histogram overflow'
   for left,right,n in zip(edges[:-1],edges[1:],counts):hist.append(dict(treatment=tr,celltype=ct,population=pop,left=left,right=right,n=int(n)))
   for pid,a in u.iterrows():donor.append(dict(treatment=tr,celltype=ct,population=pop,pid=pid,n=int(a['size']),mean=a['mean']))
 for tr in ['Untreated','CRT','CRTl','Other']:
  a=epi[(epi.treatment==tr)&(epi.level2=='Malignant')]
  for ref in ['Ductal','Ductal (atypical)']:
   b=epi[(epi.treatment==tr)&(epi.level2==ref)]
   for pop in ['All_nuclei','LYPLA1_positive_nuclei']:
    aa=a if pop=='All_nuclei' else a[a.detected];bb=b if pop=='All_nuclei' else b[b.detected]
    for unit in ['pooled_nucleus','paired_author_pid']:
     row=dict(cancer='PDAC',cohort='GSE202051',stage_id='06_EXTERNAL',run_id=RUN.name,analysis_version='gse202051_lypla1_v1',analysis_type='Malignant_vs_'+ref,
       metabolite_key='NA',metabolite_name='NA',gene='LYPLA1',unit=unit,n=len(aa),n_reference=len(bb),effect_type='rank_biserial' if unit=='pooled_nucleus' else 'mean_log1p_difference',
       effect=np.nan,ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,q_value=np.nan,test_family=unit+'__'+('untreated_primary' if tr=='Untreated' and ref=='Ductal' else 'supplementary'),
       family_n_evaluable=0,status='NOT_EVALUABLE',reason='',source_id='GEO_GSE202051_author_h5ad',treatment=tr,reference=ref,population=pop,mean_test=np.nan,mean_reference=np.nan)
     if unit=='pooled_nucleus':
      x=aa.expression.to_numpy();y=bb.expression.to_numpy()
      if min(len(x),len(y))>=20:
       u,p=stats.mannwhitneyu(x,y,alternative='two-sided',method='asymptotic');row.update(status='DONE',p_value=p,effect=2*u/(len(x)*len(y))-1,mean_test=x.mean(),mean_reference=y.mean(),reason='Exploratory; within-donor dependence unadjusted')
      else:row['reason']='Fewer than20 nuclei in a comparison group'
     else:
      x=aa.groupby('pid').expression.agg(['mean','size']);y=bb.groupby('pid').expression.agg(['mean','size']);x=x[x['size']>=20];y=y[y['size']>=20]
      ids=sorted(set(x.index)&set(y.index));x=x.loc[ids,'mean'].to_numpy();y=y.loc[ids,'mean'].to_numpy();row.update(n=len(ids),n_reference=len(ids))
      if len(ids)>=5:row.update(paired(x,y,SEED+len(rows)));row.update(status='DONE',reason='Paired author PID; sample libraries pooled within PID',mean_test=x.mean(),mean_reference=y.mean())
      else:row['reason']='Fewer than5 paired PID with >=20 selected-population nuclei per group'
     rows.append(row)
 r=pd.DataFrame(rows)
 for family,ix in r.groupby('test_family').groups.items():r.loc[ix,'q_value']=bh(r.loc[ix,'p_value']);r.loc[ix,'family_n_evaluable']=r.loc[ix,'p_value'].notna().sum()
 write(r,'results.tsv');write(pd.DataFrame(profiles),'expression_summary.tsv');write(pd.DataFrame(hist),'histogram_source.tsv')
 pd.DataFrame(donor).to_csv(PRIV/'donor_values.tsv',sep='\t',index=False)
 write(c.groupby(['treatment','level1','level2']).agg(n_cells=('cell','size'),n_donors=('pid','nunique'),n_samples=('sampleid','nunique')).reset_index(),'annotation_summary.tsv')
 print(r[r.treatment=='Untreated'][['reference','population','unit','n','n_reference','effect','p_value','q_value','mean_test','mean_reference']].to_string(index=False),flush=True)
 return r
def main():
 PUB.mkdir(exist_ok=True);PRIV.mkdir(exist_ok=True)
 spec=json.loads((RUN/'gse202051_lypla1_spec.json').read_text());(PUB/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
 if (PRIV/'cell_values.tsv.gz').exists():c=pd.read_csv(PRIV/'cell_values.tsv.gz',sep='\t');checks=json.loads((PRIV/'normalization_checks.json').read_text())
 else:c,checks=extract()
 analyze(c)
 print('HASHING_INPUT',flush=True);manifest=[dict(source_id=p.name,server_path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in [INPUT,Path(__file__),RUN/'gse202051_lypla1_spec.json']]
 write(pd.DataFrame(manifest),'source_manifest.tsv')
 (PUB/'validation.json').write_text(json.dumps(dict(status='PASS',n_nuclei=len(c),n_author_pid=c.pid.nunique(),n_libraries=c.sampleid.nunique(),
   gene_unique=True,all_count_values_nonnegative_integer=True,unique_cell_keys=True,normalization_checks=checks,author_label_consistency=True,
   python=platform.python_version(),h5py=h5py.__version__,numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,
   limitations=['Single nuclei not whole cells','Author malignancy labels reused; no new CNV validation','Nonmalignant ductal cells from tumor specimens are not healthy controls','Pooled-nucleus inference has donor dependence','No author malignant subtype assignment in this batch','Cohort independence from other studies not individually audited']),indent=2))
 (RUN/'ANALYSIS_DONE.json').write_text(json.dumps({'status':'DONE'}));print('DONE',flush=True)
if __name__=='__main__':main()
