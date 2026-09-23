"""Server-side donor-balanced disjoint metacells; not the hdWGCNA package."""
import os
os.environ['OPENBLAS_NUM_THREADS']='4'
os.environ['OMP_NUM_THREADS']='4'
from pathlib import Path
import sys,json,hashlib,importlib.util
import numpy as np,pandas as pd,h5py
from scipy import sparse
from sklearn.decomposition import PCA
R=Path(sys.argv[1]);P=R/'public';D=R/'private'
assert (R/'.running').read_text().strip()=='malignant_wgcna_v1'
P.mkdir(exist_ok=True);D.mkdir(exist_ok=True)
s=importlib.util.spec_from_file_location('helper',R/'brca_sc117_profile_v1.py');h=importlib.util.module_from_spec(s);s.loader.exec_module(h)
seed=20260923;rng=np.random.default_rng(seed)
spec=dict(version='malignant_wgcna_v1',source_commit='912db19',scope='author-labelled malignant epithelium only',method='custom donor-balanced disjoint KNN metacells followed by WGCNA;not hdWGCNA implementation',discovery='Wu2021',external='Pal2021_reprocessed',min_cells_per_donor=120,candidate_cell_cap_per_donor=1000,metacells_per_donor=8,cells_per_metacell=15,cell_overlap=0,normalization='sum raw metacell counts,all-gene library normalize10000,log1p',gene_filter='unique symbol,nonmitochondrial,detected>=5% balanced candidate cells and >=3 donors;top6000 log-normalized variance;no forced LYPLA1',pca='top2000 variable eligible Wu genes;50 PCs;study-specific PCA;within-donor neighbors',seed=seed,biological_unit='author donor labels,not independently genotype verified',inference='metacells not independent patient replicates;no cell-level P',signed_network=True,correlation='bicor,maxPOutliers=.05',soft_power_rule='smallest candidate with signed scale-free R2>=.8 and negative slope and mean connectivity>=5;else largest fit among mean connectivity>=5;fallback labelled',powers=list(range(1,11))+list(range(12,31,2)),modules=dict(minModuleSize=30,deepSplit=2,mergeCutHeight=.25,pamRespectsDendro=False),stability='leave-one-donor-out fixed-module LYPLA1 kME,not full network reclustering stability',enrichment='GO BP hypergeometric all modules;BH across all tested module-term pairs;network gene universe',new_patient_P_q=False)
(P/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
manifest=[];coverage=[];dataset={}
for co,file,cfgfile in [('Wu2021','wu_curated.h5ad','brca_sc117_Wu2021_config.json'),('Pal2021_reprocessed','pal_atlas_container.h5ad','brca_sc117_Pal2021_config.json')]:
 cfg=json.loads((R/cfgfile).read_text());src=R.parent/'20260919T140105Z_scRNA117_v1/source'/file
 manifest.append(dict(path=str(src),sha256=h.sha256(src),kind='source_h5ad'))
 with h5py.File(src,'r') as f:
  obs=h.dataframe(f['obs']);var=h.dataframe(f['raw/var']);names=var.feature_name.astype(str)
  select=obs[cfg['celltype_column']].map(cfg['celltype_map']).eq('Malignant_epithelial')
  if cfg.get('filter_column'):select &=obs[cfg['filter_column']].eq(cfg['filter_value'])
  counts=obs.loc[select,cfg['donor_column']].value_counts();keepdonors=sorted(counts[counts>=120].index.astype(str))
  sample=[]
  for donor in keepdonors:
   ix=np.flatnonzero((select&obs[cfg['donor_column']].eq(donor)).values)
   sample.extend(rng.choice(ix,min(len(ix),1000),replace=False))
  sample=np.sort(sample);select2=np.zeros(len(obs),bool);select2[sample]=True;meta=obs.iloc[sample].copy()
  raw=f['raw/X'];ptr=raw['indptr'][:];parts=[]
  for start in range(0,len(obs),2000):
   end=min(start+2000,len(obs));mask=select2[start:end]
   if mask.any():
    lo,hi=ptr[start],ptr[end];mat=sparse.csr_matrix((raw['data'][lo:hi],raw['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(var)))
    parts.append(mat[mask])
  x=sparse.vstack(parts).tocsr();assert x.shape[0]==len(meta)
  lib=np.asarray(x.sum(axis=1)).ravel();assert (lib>0).all()
  norm=x.multiply(10000/lib[:,None]).tocsr();norm.data=np.log1p(norm.data)
  mean=np.asarray(norm.mean(axis=0)).ravel();variance=np.asarray(norm.power(2).mean(axis=0)).ravel()-mean**2
  detection=np.asarray((x>0).mean(axis=0)).ravel()
  present=np.zeros(x.shape[1],int)
  for donor in keepdonors:present+=np.asarray(x[meta[cfg['donor_column']].eq(donor).values].sum(axis=0)).ravel()>0
  eligible=(~names.duplicated(keep=False)).values & ~names.str.startswith('MT-').values &(detection>=.05)&(present>=3)&(variance>0)
  if co=='Wu2021':
   gi=np.flatnonzero(eligible);gi=gi[np.argsort(-variance[gi],kind='stable')[:6000]];genes=names.iloc[gi].tolist();pcgenes=genes[:2000]
   pd.DataFrame(dict(gene=names,eligible=eligible,detection=detection,donors_detected=present,variance=variance,network_selected=names.isin(genes))).to_csv(P/'Wu_gene_filter.tsv',sep='\t',index=False)
  lookup={g:int(np.flatnonzero(names.values==g)[0]) for g in genes if (names==g).sum()==1}
  pd.DataFrame([dict(gene=g,n_source_rows=int((names==g).sum()),status='DONE' if g in lookup else 'NOT_EVALUABLE',reason='' if g in lookup else 'no unique symbol in external source') for g in genes]).to_csv(P/(co+'_gene_coverage.tsv'),sep='\t',index=False)
  available_genes=[g for g in genes if g in lookup]
  pix=[lookup[g] for g in pcgenes if g in lookup];z=norm[:,pix].toarray();z=(z-z.mean(axis=0))/np.maximum(z.std(axis=0),1e-6)
  pcs=PCA(n_components=50,svd_solver='randomized',random_state=seed).fit_transform(z)
  rows=[];mr=[];cellrows=[];donorexpr=[];donormeta=[]
  for di,donor in enumerate(keepdonors):
   available=np.flatnonzero(meta[cfg['donor_column']].eq(donor).values)
   donorraw=np.asarray(x[available].sum(axis=0)).ravel();donorexpr.append(np.log1p(donorraw[[lookup[g] for g in available_genes]]/donorraw.sum()*10000))
   sm=meta.iloc[available[0]];sub=str(sm[cfg['stratum_column']]) if cfg.get('stratum_column') else 'not_available';treat=str(sm[cfg['treatment_column']]) if cfg.get('treatment_column') else 'not_available'
   donormeta.append(dict(donor=donor,subtype=sub,treatment=treat))
   for k in range(8):
    center=rng.choice(available);dist=((pcs[available]-pcs[center])**2).sum(axis=1);chosen=available[np.argsort(dist,kind='stable')[:15]]
    available=available[~np.isin(available,chosen)]
    total=np.asarray(x[chosen].sum(axis=0)).ravel();rows.append(np.log1p(total[[lookup[g] for g in available_genes]]/total.sum()*10000));mid=f'{co}_D{di+1}_M{k+1}'
    mr.append(dict(metacell=mid,donor=donor,subtype=sub,treatment=treat,n_cells=15,library=int(total.sum())))
    cellrows.extend(dict(metacell=mid,barcode=str(meta.index[c]),donor=donor) for c in chosen)
  cellrows=pd.DataFrame(cellrows);assert not cellrows.barcode.duplicated().any()
  pd.DataFrame(rows,index=[m['metacell'] for m in mr],columns=available_genes).to_csv(D/(co+'_expr.tsv'),sep='\t',index_label='metacell')
  pd.DataFrame(mr).to_csv(D/(co+'_metacells.tsv'),sep='\t',index=False)
  cellrows.to_csv(D/(co+'_cell_membership.tsv'),sep='\t',index=False)
  pd.DataFrame(donorexpr,index=keepdonors,columns=available_genes).to_csv(D/(co+'_donor_expr.tsv'),sep='\t',index_label='donor')
  pd.DataFrame(donormeta).to_csv(D/(co+'_donor_metadata.tsv'),sep='\t',index=False)
  coverage.append(dict(cohort=co,total_malignant_cells=int(select.sum()),source_donors=len(counts),eligible_donors=len(keepdonors),candidate_cells=len(meta),network_cells=len(cellrows),metacells=len(rows),genes=len(available_genes),LYPLA1_in_network='LYPLA1' in available_genes,cell_overlap=0))
  print(coverage[-1],flush=True)
pd.DataFrame(coverage).to_csv(P/'coverage.tsv',sep='\t',index=False)
pd.DataFrame(manifest).to_csv(P/'source_manifest.tsv',sep='\t',index=False)
(P/'preparation_validation.json').write_text(json.dumps(dict(status='PASS',zero_overlap=True,equal_metacells_per_donor=8,all_raw_data_remain_server=True),indent=2))
