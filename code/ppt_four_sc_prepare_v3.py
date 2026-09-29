"""Read existing four-cancer displays and donor summaries; no inferential tests."""
from pathlib import Path
import sys,json,hashlib
import numpy as np,pandas as pd,h5py
R=Path(sys.argv[1]);P=R/'public';P.mkdir(exist_ok=True)
B=R.parents[4];C=B/'results/collaborative';sources=[];summaries=[]
def source(p,role):
 p=Path(p);sources.append(dict(path=str(p),role=role,bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest() if p.stat().st_size<100_000_000 else 'large unchanged source; original manifest retained',hash_status='DIRECT' if p.stat().st_size<100_000_000 else 'REFER_TO_PRIOR_MANIFEST'))
def read(p,role='cached values'):
 source(p,role);return pd.read_csv(p,sep='\t')
def col(g,k):
 x=g[k]
 def dec(a):return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in a])
 if isinstance(x,h5py.Group):
  cats=dec(x['categories'][:]);codes=x['codes'][:];return np.where(codes>=0,cats[np.maximum(codes,0)],'Unannotated')
 vals=x[:]
 if '__categories' in g and k in g['__categories']:
  cats=dec(g['__categories'][k][:]);return np.where(vals>=0,cats[np.maximum(vals,0)],'Unannotated')
 return dec(vals)
def makepairs(d,idcol,groupcol,value,ncol,ref,test):
 d=d[d[ncol]>=20];a=d[d[groupcol]==test].set_index(idcol)[value];b=d[d[groupcol]==ref].set_index(idcol)[value]
 assert a.index.is_unique and b.index.is_unique
 ids=sorted(set(a.index)&set(b.index));return np.c_[b.loc[ids].to_numpy(),a.loc[ids].to_numpy()]
def save(c,xy,expr,types,target,pairs,meta):
 assert len(xy)==len(expr)==len(types) and np.isfinite(xy).all() and np.isfinite(expr).all()
 delta=pairs[:,1]-pairs[:,0];meta.update(cancer=c,cells=len(xy),n_pairs=len(pairs),higher=int((delta>0).sum()),lower=int((delta<0).sum()),equal=int((delta==0).sum()),mean_difference=float(delta.mean()),target=target,cap=float(np.quantile(expr[expr>0],.99)))
 np.savez_compressed(R/(c+'_private.npz'),xy=xy,expr=expr,types=np.asarray(types,dtype=str),pairs=pairs)
 (R/(c+'_meta.json')).write_text(json.dumps(meta,ensure_ascii=False,indent=2));summaries.append(meta)
 pd.Series(types).value_counts().rename_axis('celltype').reset_index(name='cells').to_csv(P/(c+'_coverage.tsv'),sep='\t',index=False)
 print(c,json.dumps(meta,ensure_ascii=False),flush=True)
# BRCA: full-cell display cache from immediately preceding figure.
r=C/'BRCA/A/20260929T065300Z_brca_allcells_slide_v2';source(r/'private_display_cache.npz','exact preceding all-cell display')
z=np.load(r/'private_display_cache.npz');d=read(C/'BRCA/A/20260923T040012Z_epithelial_paired20_v1/private/paired_measurements.tsv')
d=d[(d.cohort=='Wu2021')&(d.partition=='ALL')&(d.gene=='LYPLA1')];pairs=d[['nonmalignant','malignant']].to_numpy()
frozen=read(r/'frozen_results.tsv','frozen donor tests');v=frozen[(frozen.cohort=='Wu2021')&(frozen.analysis_type=='ALL')&(frozen.gene=='LYPLA1')].iloc[0]
assert len(pairs)==v.n and np.isclose(np.diff(pairs).mean(),v.effect)
save('BRCA',z['xy'],z['expr'],z['types'],'Cancer Epithelial',pairs,dict(cohort='Wu 2021 · GSE176078',scope='全部细胞',p=float(v.p_value),status='FROZEN_PAIRED_TEST',ref_label='非恶性上皮',test_label='恶性上皮',locator='Cancer epithelial',note='参照为肿瘤内非恶性上皮；Naive 子集 6/7 同向，P=0.078。',conclusion='同患者比较支持恶性上皮上调。',coordinates='source X_umap; unchanged'))
frozen[(frozen.gene=='LYPLA1')&(frozen.cohort=='Wu2021')].to_csv(P/'BRCA_frozen_results.tsv',sep='\t',index=False)
# COAD: all cells plus retained CNA/CNN identity (not equating tumor tissue with malignancy).
r=C/'COAD/B';d=read(r/'20260928T031628Z_lypla1_umap_v1/private_lypla1_overlay.tsv')
e=read(r/'20260928T034349Z_lypla1_epithelial_umap_v1/private_epithelial_overlay.tsv');assert e.cell.is_unique
g=d.cell.map(e.set_index('cell').group).fillna('unclassified')
mapping=read(r/'20260928T025513Z_uhlitz_umap_v1/public/annotation_mapping.tsv').set_index('original_label').display_label
original=np.where(d.main_cell_type=='Epithelial','Epithelial',np.where(d.main_cell_type=='Immune',d.cell_type_imm_simple,d.cell_type_str_simple))
types=pd.Series(original).map(mapping).to_numpy(dtype=object);assert pd.notna(types).all()
for k,v in [('tumor_CNA','CNA epithelial'),('tumor_CNN','CNN epithelial'),('normal_reference','Normal epithelial')]:types[g==k]=v
types[(d.main_cell_type=='Epithelial')&~g.isin(['tumor_CNA','tumor_CNN','normal_reference'])]='Epithelial unclassified'
ds=read(r/'20260925T142407Z_lypla1_cell_pooled_v1/private_donor_summary.tsv');pairs=makepairs(ds,'patient','group','expression','n_cells','normal_reference','tumor_CNA')
save('COAD',d[['UMAP1','UMAP2']].to_numpy(),d.expression.to_numpy(),types,'CNA epithelial',pairs,dict(cohort='Uhlitz · GSE166555',scope='全部细胞 · 肿瘤及正常组织',p=None,status='PAIRED_DESCRIPTIVE_NO_TEST',ref_label='正常组织上皮',test_label='肿瘤 CNA 上皮',locator='CNA epithelial',note='CNA 为作者拷贝数异常标签；CNN 不等于已证实非恶性。患者图为描述性比较，未新增检验。',conclusion='CNA 上皮与正常组织上皮：同患者方向作为描述性证据。',coordinates='frozen project UMAP; no batch correction; no new embedding'))
# PDAC: all cell classes in untreated primary samples; same scope as paired evidence.
r=C/'PDAC/B/20260928T031629Z_gse202051_lypla1_v1';d=read(r/'private/cell_values.tsv.gz');d=d[d.treatment=='Untreated'].copy()
hp=B/'data/candidates/SCP1089/GSE202051_totaldata-final-toshare.h5ad';source(hp,'source author X_umap')
with h5py.File(hp) as h:
 ids=col(h['obs'],'_index');ix=pd.Index(ids).get_indexer(d.cell);assert (ix>=0).all();xy=h['obsm/X_umap'][:][ix]
types=d.level1.to_numpy(dtype=object);types[d.level2=='Malignant']='Malignant epithelial';types[d.level2=='Ductal']='Ductal';types[d.level2=='Ductal (atypical)']='Atypical ductal'
ds=read(r/'private/donor_values.tsv');ds=ds[(ds.treatment=='Untreated')&(ds.population=='All_nuclei')];pairs=makepairs(ds,'pid','celltype','mean','n','Ductal','Malignant')
frozen=read(r/'public/results.tsv','frozen donor tests');v=frozen[(frozen.unit=='paired_author_pid')&(frozen.treatment=='Untreated')&(frozen.population=='All_nuclei')&(frozen.reference=='Ductal')].iloc[0]
assert len(pairs)==v.n and np.isclose(np.diff(pairs).mean(),v.effect) and (np.diff(pairs)>0).sum()==v.higher_donors
save('PDAC',xy,d.expression.to_numpy(),types,'Malignant epithelial',pairs,dict(cohort='Hwang · GSE202051',scope='未治疗标本 · 全部细胞类型（单核）',p=float(v.p_value),status='FROZEN_PAIRED_TEST',ref_label='非恶性导管',test_label='恶性上皮',locator='Malignant epithelial',note='9 位同患者比较；全部细胞核均纳入均值。仅阳性核的既有敏感性分析不支持同向升高。',conclusion='未治疗样本的全核均值支持恶性上皮上调。',coordinates='source author X_umap; unchanged'))
frozen[(frozen.unit=='paired_author_pid')&(frozen.treatment=='Untreated')&(frozen.reference=='Ductal')].to_csv(P/'PDAC_frozen_results.tsv',sep='\t',index=False)
# PRAD: all cells in source (tumor + adjacent), paired display restricted to tumor epithelium.
r=C/'PRAD/B';d=read(r/'20260928T030756Z_lypla1_umap_v1/private/plot_values.tsv.gz').set_index('cell')
hp=r/'20260922T142000Z_discovery_v1/source/PRAD24_cellxgene.h5ad';source(hp,'source author X_umap and metadata')
with h5py.File(hp) as h:
 o=h['obs'];ids=col(o,o.attrs.get('_index','_index'));d=d.loc[ids];xy=h['obsm/X_umap'][:];assert np.allclose(xy,d[['UMAP1','UMAP2']])
 major=col(o,'celltype_major_v2');mal=col(o,'malignant_anno_merged');tissue=col(o,'type');donor=col(o,'donor_id')
types=major.astype(object)
for k,v in [('malignant','Malignant epithelial'),('normal','Normal epithelial'),('altered_benign','Altered benign epithelial')]:types[(major=='Epithelial')&(mal==k)]=v
keep=(tissue=='cancer')&(major=='Epithelial')&np.isin(mal,['normal','malignant']);ds=pd.DataFrame(dict(donor=donor[keep],group=mal[keep],expression=d.LYPLA1_log1p_cp10k.to_numpy()[keep])).groupby(['donor','group']).expression.agg(['mean','size']).reset_index()
pairs=makepairs(ds,'donor','group','mean','size','normal','malignant')
save('PRAD',xy,d.LYPLA1_log1p_cp10k.to_numpy(),types,'Malignant epithelial',pairs,dict(cohort='PRAD24 · cellxgene',scope='全部细胞 · 肿瘤及邻近组织',p=None,status='PAIRED_DESCRIPTIVE_NO_TEST',ref_label='非恶性上皮',test_label='恶性上皮',locator='Malignant epithelial',note='患者比较仅用肿瘤组织内上皮；altered benign 不纳入两组比较。描述性配对展示，未新增检验。',conclusion='保留同患者方向及反向结果；当前为描述性配对证据。',coordinates='source author X_umap; unchanged'))
pd.DataFrame(summaries).to_csv(P/'patient_summary.tsv',sep='\t',index=False)
pd.DataFrame(sources).to_csv(P/'source_manifest.tsv',sep='\t',index=False)
(P/'analysis_spec.json').write_text(json.dumps(dict(version='four_sc_v3',inference='No new tests. BRCA/PDAC existing p reused; COAD/PRAD paired descriptive summaries only.',unit='source donor labels; >=20 cells each class',color='shared viridis palette; cohort-specific positive99 cap, shown per panel; colors not compared across cohorts',data_residency='raw values and per-donor data remain server165',new_embedding=False),indent=2))
