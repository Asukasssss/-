"""Extract two study-specific single-cell inputs on server only; stream GBmap ranges."""
from pathlib import Path
import numpy as np,pandas as pd,h5py,fsspec,json
from scipy import sparse
O=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/GBM/A/20260925T103000Z_discovery_v1')
def read(q):
 if isinstance(q,h5py.Group):
  c=np.array([x.decode() if isinstance(x,bytes) else x for x in q['categories'][:]]);codes=q['codes'][:];return np.array([c[i] if i>=0 else '' for i in codes])
 return np.array([x.decode() if isinstance(x,bytes) else x for x in q[:]])
def extract(h,study,keep,group,origin):
 gp=O/'source/genes_unique.tsv';genes=pd.read_csv(gp,sep='\t');var=h[group+'var'];names=read(var['feature_name']);ensembl=read(var['_index']);idx=[];coverage=[]
 for _,g in genes.iterrows():
  loc=np.flatnonzero(names==g.gene);idx.append(int(loc[0]) if len(loc)==1 else -1);coverage.append(dict(study=study,gene=g.gene,stable_gene_id=g.stable_gene_id,source_gene_id=ensembl[loc[0]] if len(loc)==1 else '',status='DONE' if len(loc)==1 else 'NOT_EVALUABLE',reason='exact_unique_symbol' if len(loc)==1 else 'missing_or_ambiguous_symbol'))
 rows=np.flatnonzero(keep);counts=np.full((len(rows),len(genes)),np.nan,dtype=np.float32);libs=np.zeros(len(rows));good=np.array(idx)>=0;cols=np.array(idx)[good];ptr=h[group+'X/indptr'][:];shape=h[group+'X'].attrs['shape'];integer=True
 for start in range(0,len(rows),512):
  rr=rows[start:start+512];lo=int(rr.min());hi=int(rr.max())+1;b=int(ptr[lo]);e=int(ptr[hi]);vals=h[group+'X/data'][b:e];inds=h[group+'X/indices'][b:e];integer=integer and bool(np.all(vals==np.floor(vals)));assert np.all(vals>=0)
  mat=sparse.csr_matrix((vals,inds,ptr[lo:hi+1]-b),shape=(hi-lo,int(shape[1])))[rr-lo];libs[start:start+len(rr)]=np.asarray(mat.sum(axis=1)).ravel();counts[start:start+len(rr),good]=mat[:,cols].toarray()
  if start%4096==0:print(study,'cells_read',start+len(rr),flush=True)
 assert np.all(libs>0);z=np.log1p(counts*10000/libs[:,None]).astype(np.float32);obs=pd.DataFrame({k:read(h['obs/'+k])[rows] for k in ['donor_id','cell_type']});obs['location']=read(h['obs/Location' if study=='Darmanis2017' else 'obs/location'])[rows];obs['annotation_origin']=origin
 np.savez_compressed(O/'private'/f'{study}_candidate_expression.npz',z=z,detected=counts>0,measured=good,genes=genes.gene.to_numpy(str),library=libs)
 obs.to_csv(O/'private'/f'{study}_cell_metadata.tsv',sep='\t',index=False);pd.DataFrame(coverage).to_csv(O/'source'/f'{study}_gene_coverage.tsv',sep='\t',index=False)
 summary=dict(study=study,cells=len(obs),source_donor_labels=obs.donor_id.nunique(),gene_library_size=int(shape[1]),genes_measured=int(good.sum()),candidate_genes=len(genes),nonnegative_integer_input=integer,normalization='log1p(10000*gene_count/full_retained_gene_library);Smartseq_read_counts_not_UMI',annotation_origin=origin,source_group=group or 'root',source_tissue_scope='GBM tumors plus original sampled margins;not healthy donor reference')
 (O/'source'/f'{study}_extraction_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
if not (O/'private/Darmanis2017_candidate_expression.npz').exists():
 with h5py.File(O/'source/Darmanis2017_cellxgene.h5ad','r') as h:extract(h,'Darmanis2017',np.ones(len(h['obs/_index']),bool),'','author_clusters_standardized_by_CELLxGENE')
if not (O/'private/Neftel2019_candidate_expression.npz').exists():
 with fsspec.open('https://datasets.cellxgene.cziscience.com/861acfd8-25f0-418b-a445-aa96da232827.h5ad','rb',block_size=4*2**20).open() as f:
  with h5py.File(f,'r') as h:extract(h,'Neftel2019',read(h['obs/author'])=='Neftel2019','raw/','GBmap_harmonized_reannotation_of_Neftel2019')
