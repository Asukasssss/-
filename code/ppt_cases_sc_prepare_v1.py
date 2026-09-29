"""Server-only extraction for ASNS, UCKL1, SLC6A6; no new tests."""
from pathlib import Path
import sys,json
import numpy as np,pandas as pd,h5py
from scipy import sparse
R=Path(sys.argv[1]);P=R/'public';P.mkdir(exist_ok=True)
B=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716');C=B/'results/collaborative';OLD=C/'BRCA/A/20260929T071500Z_four_sc_v3'
def dec(a):return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in a])
def col(g,k):
 x=g[k]
 if isinstance(x,h5py.Group):
  cats=dec(x['categories'][:]);codes=x['codes'][:];return np.where(codes>=0,cats[np.maximum(codes,0)],'Unannotated')
 a=x[:]
 if '__categories' in g and k in g['__categories']:return dec(g['__categories'][k][:])[a]
 return dec(a)
def expression(path,gene,counts='raw/X',var='raw/var',name='feature_name'):
 cache=R/(gene+'_'+path.name+'.npz')
 if cache.exists():z=np.load(cache);return z['ids'],z['expr']
 with h5py.File(path) as h:
  ids=col(h['obs'],h['obs'].attrs.get('_index','_index'));names=col(h[var],name);ix=np.flatnonzero(names==gene);assert len(ix)==1
  a=h[counts];ptr=a['indptr'][:];expr=np.empty(len(ids))
  for start in range(0,len(ids),4000):
   end=min(start+4000,len(ids));lo,hi=ptr[start],ptr[end]
   x=sparse.csr_matrix((a['data'][lo:hi],a['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,len(names)))
   assert np.isfinite(x.data).all() and (x.data>=0).all() and np.allclose(x.data,np.round(x.data))
   total=np.asarray(x.sum(1)).ravel();assert (total>0).all();expr[start:end]=np.log1p(x[:,ix[0]].toarray().ravel()/total*10000)
 np.savez_compressed(cache,ids=ids,expr=expr);print('EXTRACTED',gene,len(ids),flush=True);return ids,expr
def pairs(df,ref,test):
 a=df.groupby(['donor','group']).expression.agg(['mean','size']).reset_index();a=a[a['size']>=20]
 n=a[a.group==ref].set_index('donor')['mean'];t=a[a.group==test].set_index('donor')['mean'];common=sorted(set(n.index)&set(t.index));return np.c_[n.loc[common],t.loc[common]]
summaries=[]
def save(c,g,expr,ab,meta_updates):
 z=np.load(OLD/(c+'_private.npz'));m=json.loads((OLD/(c+'_meta.json')).read_text());assert len(expr)==len(z['expr'])
 m.update(meta_updates);d=ab[:,1]-ab[:,0];m.update(gene=g,n_pairs=len(ab),higher=int((d>0).sum()),lower=int((d<0).sum()),equal=int((d==0).sum()),mean_difference=float(d.mean()),cap=float(np.quantile(expr[expr>0],.99)))
 np.savez_compressed(R/(c+'_private.npz'),xy=z['xy'],types=z['types'],expr=expr,pairs=ab)
 (R/(c+'_meta.json')).write_text(json.dumps(m,ensure_ascii=False,indent=2));summaries.append(m);print('PREPARED',c,g,'pairs',len(ab),'up',m['higher'],flush=True)
 for k in ['cells','cancer']:assert k in m
 pd.Series(z['types']).value_counts().rename_axis('celltype').reset_index(name='cells').to_csv(P/(c+'_coverage.tsv'),sep='\t',index=False)
 pd.DataFrame({'type':z['types'],'expression':expr}).groupby('type').expression.agg(['size','mean']).to_csv(P/(c+'_'+g+'_celltype_descriptive.tsv'),sep='\t')
# ASNS: exact Wu all-cell coordinates, paired comparison derived descriptively.
hp=C/'BRCA/A/20260919T140105Z_scRNA117_v1/source/wu_curated.h5ad';ids,x=expression(hp,'ASNS')
with h5py.File(hp) as h:donor=col(h['obs'],'donor_id');typ=col(h['obs'],'celltype_major')
z=np.load(OLD/'BRCA_private.npz');assert np.array_equal(typ,z['types'])
ab=pairs(pd.DataFrame(dict(donor=donor,group=typ,expression=x)),'Normal Epithelial','Cancer Epithelial')
save('BRCA','ASNS',x,ab,dict(p=None,status='PAIRED_DESCRIPTIVE_NO_NEW_TEST',conclusion='同患者方向为描述性结果；参照为肿瘤内非恶性上皮。',note='同患者两类上皮均至少 20 个细胞；患者平均表达配对，未新增显著性检验。'))
# UCKL1: cached exact full-library raw counts and frozen CNA reference rules.
r=C/'COAD/B/20260922T141755Z_source_contract_sc_v2';z=np.load(r/'private_cell_counts.npz');m=pd.read_csv(r/'private_cell_metadata.tsv',sep='\t');x=np.log1p(z['UCKL1']/z['total']*10000)
d=pd.read_csv(C/'COAD/B/20260928T031628Z_lypla1_umap_v1/private_lypla1_overlay.tsv',sep='\t');ix=pd.Index(m.cell_id).get_indexer(d.cell);assert (ix>=0).all();x=x[ix]
e=pd.read_csv(C/'COAD/B/20260928T034349Z_lypla1_epithelial_umap_v1/private_epithelial_overlay.tsv',sep='\t');groups=d.cell.map(e.set_index('cell').group).fillna('other')
ab=pairs(pd.DataFrame(dict(donor=m.set_index('cell_id').loc[d.cell,'patient'].to_numpy(),group=groups,expression=x)),'normal_reference','tumor_CNA')
frozen=pd.read_csv(C/'COAD/B/20260925T140257Z_uckl1_cna_v1/public/results.tsv',sep='\t');v=frozen[frozen.contrast=='CNA_vs_normal_reference'].iloc[0]
assert len(ab)==v.n and np.isclose(np.diff(ab).mean(),v.mean_difference)
save('COAD','UCKL1',x,ab,dict(p=float(v.p_value),status='FROZEN_PAIRED_TEST',conclusion='多数患者 CNA 上皮更高，但配对检验未显著。',note='CNA 标签沿用作者；正常组织上皮参照；P 复用既有精确符号检验。'))
# Hwang counts layer, all untreated nuclei exactly as LYPLA1 display.
hp=B/'data/candidates/SCP1089/GSE202051_totaldata-final-toshare.h5ad';ids,x=expression(hp,'SLC6A6','layers/counts','var','_index')
d=pd.read_csv(C/'PDAC/B/20260928T031629Z_gse202051_lypla1_v1/private/cell_values.tsv.gz',sep='\t');d=d[d.treatment=='Untreated'];ix=pd.Index(ids).get_indexer(d.cell);assert (ix>=0).all();x=x[ix]
ab=pairs(pd.DataFrame(dict(donor=d.pid.to_numpy(),group=d.level2.to_numpy(),expression=x)),'Ductal','Malignant')
save('PDAC','SLC6A6',x,ab,dict(p=None,status='PAIRED_DESCRIPTIVE_NO_NEW_TEST',conclusion='以全细胞背景定位 SLC6A6；上皮比较按同患者描述。',note='Hwang 未治疗样本；仅描述性患者比较，无新增检验。'))
# PRAD: right-hand patient comparison uses existing matched myeloid results.
hp=C/'PRAD/B/20260922T142000Z_discovery_v1/source/PRAD24_cellxgene.h5ad';ids,x=expression(hp,'SLC6A6')
d=pd.read_csv(C/'PRAD/B/20260928T030756Z_lypla1_umap_v1/private/plot_values.tsv.gz',sep='\t');assert np.array_equal(ids,d.cell.to_numpy())
r=C/'PRAD/B/20260925T133653Z_slc6a6_sc_paired_v1';abdf=pd.read_csv(r/'private/paired_donor_values_private.tsv',sep='\t');abdf=abdf[abdf.contrast=='Myeloid'];ab=abdf[['mean_log1p_cp10k_normal','mean_log1p_cp10k_tumor']].to_numpy()
v=pd.read_csv(r/'public/06_EXTERNAL/results.tsv',sep='\t');v=v[v.contrast=='Myeloid'].iloc[0];assert len(ab)==v.n and np.isclose(np.diff(ab).mean(),v.effect)
save('PRAD','SLC6A6',x,ab,dict(p=float(v.p_value),status='FROZEN_PAIRED_TEST',target='Myeloid',locator='Myeloid',ref_label='癌旁组织髓系',test_label='肿瘤组织髓系',conclusion='髓系肿瘤—癌旁配对未显著；组织 RNA 下降不能直接归于髓系。',note='同类髓系按作者供者配对；不进行恶性髓系推断。'))
pd.DataFrame(summaries).to_csv(P/'sc_patient_summary.tsv',sep='\t',index=False)
print('SC_READY',flush=True)
