"""Original-author CRA001160 LYPLA1 analysis. Execute only in a new server run."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import hashlib,json,tarfile,shutil,sys,platform
import numpy as np
import pandas as pd
import scipy
from scipy import io,stats
RUN=Path(__file__).resolve().parent
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
DATA=ROOT/'data/candidates/CRA001160_author_v1'
PUB=RUN/'public';PRIVATE=RUN/'private'
GENE='LYPLA1';SEED=20260927
D1='Ductal cell type 1';D2='Ductal cell type 2'
CONTRASTS=['Within_tumor_type2_vs_type1_paired','Tumor_type2_vs_control_type1_unpaired']
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def write(df,name):df.to_csv(PUB/name,sep='\t',index=False,na_rep='NA')
def bh(p):
 p=np.asarray(p,float);out=np.full(len(p),np.nan);ok=np.where(np.isfinite(p))[0]
 order=ok[np.argsort(p[ok])];n=len(order)
 if n:out[order]=np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1].clip(0,1)
 return out
def inference(x,y,paired,seed=SEED):
 x=np.asarray(x,float);y=np.asarray(y,float)
 if min(len(x),len(y))<5:return {'status':'NOT_EVALUABLE','reason':'Fewer than five eligible units in at least one group'}
 assert np.isfinite(x).all() and np.isfinite(y).all()
 if paired:assert len(x)==len(y)
 rng=np.random.default_rng(seed);effect=x.mean()-y.mean();extreme=0
 for _ in range(100):
  if paired:
   perm=(rng.choice([-1.,1.],size=(1000,len(x)))*(x-y)).mean(axis=1)
  else:
   z=np.r_[x,y];order=np.argsort(rng.random((1000,len(z))),axis=1)
   perm=z[order[:,:len(x)]].mean(axis=1)-z[order[:,len(x):]].mean(axis=1)
  extreme+=np.count_nonzero(np.abs(perm)>=abs(effect)-1e-12)
 if paired:
  bs=(x-y)[rng.integers(0,len(x),(20000,len(x)))].mean(axis=1)
  rankp=1. if np.all(x==y) else stats.wilcoxon(x,y,zero_method='wilcox',alternative='two-sided').pvalue
 else:
  bs=x[rng.integers(0,len(x),(20000,len(x)))].mean(axis=1)-y[rng.integers(0,len(y),(20000,len(y)))].mean(axis=1)
  rankp=stats.mannwhitneyu(x,y,alternative='two-sided',method='auto').pvalue
 lo,hi=np.quantile(bs,[.025,.975])
 return dict(status='DONE',reason='',effect=effect,ci_lower=lo,ci_upper=hi,p_value=(extreme+1)/100001,
             rank_test_p=rankp,mean_test=x.mean(),mean_reference=y.mean(),median_test=np.median(x),median_reference=np.median(y),
             n_higher=int(np.sum(x>y)) if paired else np.nan,n_lower=int(np.sum(x<y)) if paired else np.nan,
             n_equal=int(np.sum(x==y)) if paired else np.nan)
def extract():
 manifest=json.loads((DATA/'download_manifest.json').read_text())
 assert sha(DATA/'PDAC.tar.gz')==manifest['sha256']
 dest=DATA/'extracted';dest.mkdir(exist_ok=True)
 # Never trust archive paths or allow symlinks.
 with tarfile.open(DATA/'PDAC.tar.gz','r|gz') as tar:
  for item in tar:
   target=(dest/item.name).resolve()
   assert dest.resolve() in target.parents or target==dest.resolve()
   if item.isdir():target.mkdir(exist_ok=True,parents=True)
   elif item.isfile():
    if target.exists() and target.stat().st_size==item.size:continue
    target.parent.mkdir(exist_ok=True,parents=True)
    with tar.extractfile(item) as src,target.open('wb') as out:shutil.copyfileobj(src,out)
   else:raise ValueError('Unexpected archive member')
 return dest
def load_cells():
 dest=extract();ann=pd.read_csv(DATA/'all_celltype.txt',sep='\t')
 assert list(ann.columns)==['cell.name','cluster'] and ann['cell.name'].is_unique and ann.notna().all().all()
 assert len(ann)==57530
 ann['unit']=ann['cell.name'].str.split('_').str[0]
 ann=ann.set_index('cell.name');pieces=[];audit=[];seen=set()
 for folder in sorted(p for p in dest.iterdir() if p.is_dir()):
  unit=folder.name
  genes=pd.read_csv(folder/'genes.tsv',sep='\t',header=None,dtype=str)
  barcodes=pd.read_csv(folder/'barcodes.tsv',sep='\t',header=None,dtype=str).iloc[:,0]
  assert barcodes.is_unique
  gene_col=1 if genes.shape[1]>1 else 0
  hits=np.where(genes.iloc[:,gene_col].values==GENE)[0]
  assert len(hits)==1,(unit,'LYPLA1 gene label not unique')
  # GSA 10x barcodes carry -1; original annotation uses the explicit library prefix.
  cell_ids=unit+'_'+barcodes.str.replace(r'-1$','',regex=True)
  assert cell_ids.is_unique
  expected=ann.index[ann.unit==unit]
  assert set(expected).issubset(set(cell_ids)),(unit,'Missing author-annotated barcode')
  m=io.mmread(str(folder/'matrix.mtx')).tocsr()
  assert m.shape==(len(genes),len(barcodes))
  assert np.isfinite(m.data).all() and (m.data>=0).all() and np.equal(m.data,np.floor(m.data)).all()
  lib=np.asarray(m.sum(axis=0)).ravel();counts=m[hits[0],:].toarray().ravel()
  ng=np.asarray((m>0).sum(axis=0)).ravel()
  keep=cell_ids.isin(expected).values
  assert (lib[keep]>0).all()
  ids=cell_ids[keep].values;labels=ann.loc[ids,'cluster'].values
  seen.update(ids)
  df=pd.DataFrame({'cell':ids,'unit':unit,'tissue':'Tumor' if unit.startswith('T') else 'Control',
                   'cell_type':labels,'count':counts[keep],'library':lib[keep],'n_genes':ng[keep]})
  df['log1p_10k']=np.log1p(10000*df['count']/df.library);df['detected']=(df['count']>0).astype(float)
  pieces.append(df)
  audit.append(dict(unit=unit,matrix_cells=len(barcodes),author_cells=len(df),n_genes=len(genes),
    library_median=float(np.median(lib[keep])),matrix_sha256=sha(folder/'matrix.mtx'),genes_sha256=sha(folder/'genes.tsv'),barcodes_sha256=sha(folder/'barcodes.tsv')))
  print('EXTRACTED',unit,'annotated_cells',len(df),flush=True)
  del m
 assert seen==set(ann.index) and len(pieces)==35
 cells=pd.concat(pieces,ignore_index=True)
 assert cells.cell.is_unique and len(cells)==57530
 cells.to_csv(PRIVATE/'cell_values.tsv.gz',sep='\t',index=False)
 pd.DataFrame(audit).to_csv(PRIVATE/'sample_manifest.tsv',sep='\t',index=False)
 return cells
def unit_table(cells):
 u=cells.groupby(['unit','tissue','cell_type']).agg(n_cells=('count','size'),log1p_10k=('log1p_10k','mean'),
  detection=('detected','mean'),sum_gene=('count','sum'),sum_library=('library','sum'),median_library=('library','median')).reset_index()
 u['pseudobulk_log2cpm']=np.log2((u.sum_gene+.5)*1e6/u.sum_library)
 return u
def analyze(units):
 rows=[];cache={}
 for cutoff,metric in [(20,'log1p_10k'),(10,'log1p_10k'),(50,'log1p_10k'),(20,'detection'),(20,'pseudobulk_log2cpm')]:
  good=units[units.n_cells>=cutoff]
  for ci,contrast in enumerate(CONTRASTS):
   a=good[(good.tissue=='Tumor')&(good.cell_type==D2)].set_index('unit')
   b=good[(good.tissue==('Tumor' if ci==0 else 'Control'))&(good.cell_type==D1)].set_index('unit')
   paired=ci==0
   if paired:
    shared=sorted(set(a.index)&set(b.index));a=a.loc[shared];b=b.loc[shared]
   key=(ci,metric,tuple(a.index),tuple(b.index),tuple(a[metric]),tuple(b[metric]))
   reused=key in cache
   if not reused:cache[key]=inference(a[metric].values,b[metric].values,paired,SEED+cutoff+ci)
   res=dict(cache[key])
   row=dict(cancer='PDAC',cohort='CRA001160',stage_id='06_EXTERNAL',run_id=RUN.name,analysis_version='cra001160_lypla1_v1',
    analysis_type=contrast,metabolite_key='NA',metabolite_name='NA',gene=GENE,unit='author_patient_sample_label',n=len(a),n_reference=len(b),
    effect_type=metric+'_mean_difference',effect=np.nan,ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,q_value=np.nan,
    test_family='LYPLA1_two_primary_contrasts' if cutoff==20 and metric=='log1p_10k' else 'sensitivity_separate_descriptive',
    family_n_evaluable=np.nan,status=res['status'],reason=res['reason'],source_id='GSA_CRA001160_original_processed',
    min_cells=cutoff,metric=metric,cells_test=int(a.n_cells.sum()),cells_reference=int(b.n_cells.sum()),paired=paired,
    identical_input_reused=reused)
   row.update(res);rows.append(row)
 out=pd.DataFrame(rows);primary=out.test_family=='LYPLA1_two_primary_contrasts'
 out.loc[primary,'q_value']=bh(out.loc[primary,'p_value'])
 out.loc[primary,'family_n_evaluable']=out.loc[primary,'p_value'].notna().sum()
 write(out,'results.tsv');return out
def summarize(cells,units):
 rng=np.random.default_rng(SEED);rows=[]
 for (tissue,ct),g in units.groupby(['tissue','cell_type']):
  eligible=g[g.n_cells>=20];d=dict(tissue=tissue,cell_type=ct,n_cells=int(g.n_cells.sum()),n_units_any=len(g),
   n_units_eligible=len(eligible),n_cells_eligible=int(eligible.n_cells.sum()))
  for metric in ['log1p_10k','detection','median_library']:
   x=eligible[metric].values
   if len(x):
    bs=x[rng.integers(0,len(x),(20000,len(x)))].mean(axis=1);lo,hi=np.quantile(bs,[.025,.975])
    d.update({metric+'_mean':float(x.mean()),metric+'_ci_lower':lo,metric+'_ci_upper':hi})
  rows.append(d)
 write(pd.DataFrame(rows),'celltype_summary.tsv')
 write(cells.groupby(['tissue','cell_type']).size().rename('n_cells').reset_index(),'annotation_counts.tsv')
def main():
 PUB.mkdir(exist_ok=True);PRIVATE.mkdir(exist_ok=True)
 spec=json.loads((RUN/'cra001160_lypla1_spec.json').read_text());shutil.copy2(RUN/'cra001160_lypla1_spec.json',PUB/'analysis_spec.json')
 cells=pd.read_csv(PRIVATE/'cell_values.tsv.gz',sep='\t') if '--resume-summary' in sys.argv else load_cells()
 units=unit_table(cells);units.to_csv(PRIVATE/'unit_values.tsv',sep='\t',index=False)
 results=analyze(units);summarize(cells,units)
 manifests=[]
 for p,url in [(DATA/'PDAC.tar.gz','https://download.cncb.ac.cn/gsa/CRA001160/other_files/PDAC.tar.gz'),(DATA/'all_celltype.txt','https://download.cncb.ac.cn/gsa/CRA001160/other_files/all_celltype.txt'),(RUN/'cra001160_lypla1.py','code/pdac/cra001160_lypla1.py'),(RUN/'cra001160_lypla1_spec.json','code/pdac/cra001160_lypla1_spec.json')]:
  manifests.append(dict(source_id=p.name,path_or_url=url,server_path=str(p),sha256=sha(p),bytes=p.stat().st_size))
 write(pd.DataFrame(manifests),'source_manifest.tsv')
 validation=dict(status='PASS',original_author_cells=57530,matched_cells=len(cells),unique_cells=bool(cells.cell.is_unique),
  n_author_units=int(cells.unit.nunique()),n_tumor_units=int(cells[cells.tissue=='Tumor'].unit.nunique()),n_control_units=int(cells[cells.tissue=='Control'].unit.nunique()),
  nonnegative_integer_counts=True,lypla1_unique_each_matrix=True,all_gene_library_sizes=True,range_etag_and_length_verified=True,
  archive_sha256_verified=True,primary_family_size=2,per_cell_statistics_exported=False,
  python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,
  limitations=['Original author type1/type2 labels reused; no new CNV confirmation','Clinical independence and control diagnoses require original source context','No tumor-control pairing inferred from matching numeric suffixes','Cross-study independence from existing cohorts not independently established'])
 (PUB/'validation.json').write_text(json.dumps(validation,indent=2))
 (RUN/'ANALYSIS_DONE.json').write_text(json.dumps({'status':'DONE','public':str(PUB)},indent=2))
 print(results[['analysis_type','metric','min_cells','n','n_reference','effect','p_value','q_value','status']].to_string(index=False),flush=True)
if __name__=='__main__':main()
