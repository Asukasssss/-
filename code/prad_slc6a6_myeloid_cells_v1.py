"""User-requested pooled-cell SLC6A6 comparison using depositor myeloid labels.

No donor weighting/pairing. Cell-level P values are not patient-level inference.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import argparse, json, hashlib, platform
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scipy
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--source-run',type=Path,required=True);ap.add_argument('--code-commit',required=True);ap.add_argument('--author-commit',required=True);a=ap.parse_args()
out=a.out
with (out/'.running').open('x') as f:f.write('prad_slc6a6_myeloid_cells_v1')
pub=out/'public/06_EXTERNAL';pub.mkdir(parents=True,exist_ok=False)
priv=out/'private';priv.mkdir()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA')
def js(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
src=a.source_run/'source/PRAD24_cellxgene.h5ad';lp=a.source_run/'private/sc_library_private.tsv'
adata=ad.read_h5ad(src,backed='r');obs=adata.obs.copy()
assert obs.index.is_unique and obs.donor_id.notna().all()
assert obs.groupby('sample',observed=True).donor_id.nunique().max()==1
assert obs.groupby('sample',observed=True)['type'].nunique().max()==1
mask=obs.celltype_major_v2.astype(str).eq('Myeloid').to_numpy();my=obs.loc[mask].copy()
expected={'Macrophages','Monocytes','cDC2','Mast_cells','cDC1','pDC','mDC_regs'}
assert set(my.celltype_minor_v2.astype(str))==expected
assert not my[['celltype_minor_v2','celltype_subset_v2']].isna().any().any()
assert my.groupby('celltype_subset_v2',observed=True).celltype_minor_v2.nunique().max()==1
expected_sub=['Monocytes_S100A9+/c0','Macrophages_SELENOP+/c1','cDC2_CLEC10A+/c2','Macrophages_C3+/c3','cDC2_CD207+/c4','Macrophages_C1QB+/c5','Monocytes_S100A9+/c6','Monocytes_FCGR3A+/c7','Mast_cells_TPSB2+/c8','Macrophages_MT1G+/c9','cDC1_CLEC9A+/c10','pDC_IRF4+/c11','mDC_regs_BIRC3+/c12']
assert set(my.celltype_subset_v2.astype(str))==set(expected_sub)
fields=['celltype_major_v2','celltype_minor_v2','celltype_subset_v2','cell_type','cell_type_ontology_term_id']
ann=my.groupby(fields,observed=True).size().reset_index(name='n_cells')
ann['origin']='depositor v2;matched original annotation plus July2023 update;no reannotation'
save(ann,pub/'author_annotation_map.tsv')
genes=['SLC6A6','S100A9','SELENOP','CLEC10A','C3','CD207','C1QB','FCGR3A','TPSB2','MT1G','CLEC9A','IRF4','BIRC3']
v=adata.raw.var;symbols=v.feature_name.astype(str).to_numpy();ix=[]
for g in genes:
 hits=np.flatnonzero(symbols==g);assert len(hits)==1,(g,len(hits));ix.append(hits[0])
lib=pd.read_csv(lp,sep='\t',index_col=0);assert lib.index.is_unique and lib.index.equals(obs.index)
all_lib=lib.library.to_numpy(float)
norm=np.zeros((len(obs),len(genes)),np.float32);detect=np.zeros_like(norm,dtype=bool)
err=0.
for start in range(0,len(obs),2000):
 end=min(start+2000,len(obs));sel=mask[start:end]
 if not sel.any():continue
 x=adata.raw.X[start:end].tocsr()[sel];L=np.asarray(x.sum(1)).ravel()
 assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.round(x.data)).all()
 assert np.array_equal(L,all_lib[start:end][sel]) and (L>0).all()
 raw=x[:,ix].toarray();positions=np.arange(start,end)[sel]
 norm[positions]=np.log1p(10000*raw/L[:,None]);detect[positions]=raw>0
print('raw myeloid counts and full libraries verified',int(mask.sum()),flush=True)
my['SLC6A6_log1p_CP10K']=norm[mask,0];my['SLC6A6_detected']=detect[mask,0]
save(my[['donor_id','sample','type','celltype_minor_v2','celltype_subset_v2','SLC6A6_log1p_CP10K','SLC6A6_detected']],priv/'cell_values_private.tsv')
# Reproduce prior donor-level Myeloid means to detect index/normalization drift.
old=pd.read_csv(a.source_run/'private/sc_donor_profiles_private.tsv',sep='\t')
old=old[old.gene.eq('SLC6A6')&old.celltype.eq('Myeloid')]
fresh=my.groupby(['type','donor_id'],observed=True).SLC6A6_log1p_CP10K.mean().reset_index()
check=fresh.merge(old,left_on=['type','donor_id'],right_on=['partition','donor'],validate='one_to_one')
assert len(check)==len(fresh)==len(old)
maxerror=float(abs(check.SLC6A6_log1p_CP10K-check.mean_log1p_cp10k).max());assert maxerror<1e-6
rows=[]
for field,level in [('celltype_minor_v2','minor'),('celltype_subset_v2','subset')]:
 for label,z in my.groupby(field,observed=True):
  t=z[z['type'].astype(str).eq('cancer')];n=z[z['type'].astype(str).eq('adj_benign')]
  x=t.SLC6A6_log1p_CP10K.to_numpy(float);y=n.SLC6A6_log1p_CP10K.to_numpy(float)
  row=dict(cancer='PRAD',cohort='PRAD24_CELLxGENE',stage_id='06_EXTERNAL',run_id=out.name,analysis_version='prad_slc6a6_myeloid_cells_v1',analysis_type='pooled_cell_tissue_comparison',metabolite_key='NA',metabolite_name='NA',gene='SLC6A6',unit='cell_unpaired_no_donor_adjustment',n=len(x),n_reference=len(y),effect_type='Cliffs_delta_from_Mann_Whitney_U',effect=np.nan,ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,q_value=np.nan,test_family='SLC6A6_POOLED_'+level.upper(),family_n_evaluable=0,status='NOT_EVALUABLE',reason='fewer_than3_cells_in_one_tissue',source_id='PRAD24_CELLxGENE',annotation_level=level,celltype=str(label),mean_tumor=float(x.mean()) if len(x) else np.nan,mean_adjacent=float(y.mean()) if len(y) else np.nan,median_tumor=float(np.median(x)) if len(x) else np.nan,median_adjacent=float(np.median(y)) if len(y) else np.nan,detected_fraction_tumor=float(t.SLC6A6_detected.mean()) if len(x) else np.nan,detected_fraction_adjacent=float(n.SLC6A6_detected.mean()) if len(y) else np.nan,n_donors_tumor=int(t.donor_id.nunique()),n_donors_adjacent=int(n.donor_id.nunique()),inference_scope='cell-level only;cells within donor not independent')
  row['mean_difference']=row['mean_tumor']-row['mean_adjacent']
  if min(len(x),len(y))>=3:
   u,pv=stats.mannwhitneyu(x,y,alternative='two-sided',method='asymptotic',use_continuity=True)
   ranks=stats.rankdata(np.r_[x,y]);u2=float(ranks[:len(x)].sum()-len(x)*(len(x)+1)/2)
   assert abs(u-u2)<1e-8
   row.update(effect=float(2*u/(len(x)*len(y))-1),p_value=float(pv),status='DONE',reason='cell_level_exploratory',U=float(u))
  rows.append(row)
d=pd.DataFrame(rows)
for f,idx in d.groupby('test_family').groups.items():
 ok=d.index.isin(idx)&d.p_value.notna();d.loc[idx,'family_n_evaluable']=int(ok.sum())
 if ok.any():d.loc[ok,'q_value']=multipletests(d.loc[ok,'p_value'],method='fdr_bh')[1]
save(d,pub/'results.tsv')
# Marker expression is only a check of depositor labels, never used to assign new types.
markerrows=[]
subobs=obs.loc[mask].copy().reset_index(drop=True);vnorm=norm[mask];vdet=detect[mask]
for subset,idx in subobs.groupby('celltype_subset_v2',observed=True).indices.items():
 for j,g in enumerate(genes[1:],start=1):markerrows.append(dict(author_subset=str(subset),gene=g,n_cells=len(idx),mean_log1p_CP10K=float(vnorm[idx,j].mean()),detection_fraction=float(vdet[idx,j].mean()),scope='pooled descriptor;not independent subtype validation'))
save(pd.DataFrame(markerrows),pub/'author_label_marker_profiles.tsv')
markers=pd.DataFrame(markerrows).pivot(index='author_subset',columns='gene',values='mean_log1p_CP10K').reindex(expected_sub)
scaled=markers/(markers.max(axis=0).replace(0,1))
fig,ax=plt.subplots(figsize=(10,7));im=ax.imshow(scaled,aspect='auto',cmap='Blues',vmin=0,vmax=1);ax.set_yticks(range(len(scaled)));ax.set_yticklabels(scaled.index,fontsize=8);ax.set_xticks(range(len(scaled.columns)));ax.set_xticklabels(scaled.columns,rotation=60,ha='right');ax.set_title('Author myeloid labels: marker expression check\nWithin-gene scaling for display; no reannotation');fig.colorbar(im,ax=ax,label='Mean / maximum across author subsets');fig.tight_layout();fig.savefig(pub/'author_marker_check.png',dpi=180);plt.close(fig)
major=d[d.annotation_level.eq('minor')].set_index('celltype').reindex(['Monocytes','Macrophages','cDC2','cDC1','pDC','mDC_regs','Mast_cells'])
fig,axs=plt.subplots(1,2,figsize=(13,5),sharey=True);yy=np.arange(len(major));h=.35
for ax,cols,title in [(axs[0],['mean_tumor','mean_adjacent'],'Mean log1p(CP10K)'),(axs[1],['detected_fraction_tumor','detected_fraction_adjacent'],'SLC6A6 detected cells (%)')]:
 scale=100 if ax is axs[1] else 1
 ax.barh(yy-h/2,major[cols[0]]*scale,h,color='#cb6557',label='Cancer');ax.barh(yy+h/2,major[cols[1]]*scale,h,color='#438aaa',label='Adjacent benign');ax.set_xlabel(title);ax.spines[['top','right']].set_visible(False)
axs[0].set_yticks(yy);axs[0].set_yticklabels([f'{label} (T={int(r.n)}, N={int(r.n_reference)})\ncell-level P={r.p_value:.3g}' for label,r in major.iterrows()],fontsize=9);axs[0].invert_yaxis();axs[1].legend();fig.suptitle('SLC6A6 in author myeloid subtypes | pooled cells\nUnpaired cell-level comparison; no donor adjustment',fontsize=12);fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(pub/f'SLC6A6_myeloid_pooled.{ext}',dpi=180)
plt.close(fig)
spec=dict(version='prad_slc6a6_myeloid_cells_v1',parent_commit=a.code_commit,author_repository='https://github.com/swarbricklab/apostolov_pca_atlas',author_commit=a.author_commit,primary_labels='celltype_minor_v2;7 author types',secondary_labels='celltype_subset_v2;13 author subsets',selection='all observed myeloid minor/subset labels;no outcome filtering',normalization='log1p(10000*raw_count/full_raw_gene_library);same as prior',test='two-sided Mann-Whitney U;asymptotic tie and continuity corrected;zeros retained',minimum_cells_each_tissue=3,unit='cell;no donor pairing or adjustment per user request',P_threshold=.05,BH='separate seven-minor and thirteen-subset families;reported not selection gate',scientific_scope='pooled expression difference;not patient-level evidence;possible donor/cell-composition confounding',CI='NA;no cell bootstrap CI implying independent patients',software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,anndata=ad.__version__,scipy=scipy.__version__))
js(pub/'analysis_spec.json',spec)
js(pub/'validation.json',dict(status='DONE',myeloid_cells=len(my),cancer_cells=int(my['type'].eq('cancer').sum()),adjacent_cells=int(my['type'].eq('adj_benign').sum()),seven_author_minor_types=True,thirteen_author_subsets=True,subset_to_minor_unambiguous=True,exact_gene_matches=True,raw_integer_counts_and_cached_library_verified=True,prior_myeloid_profile_max_error=maxerror,independent_U_statistic_verified=True,new_clustering=False,new_annotation=False,genotype_identity_validation='NOT_RUN;source matrices only;inherit depositor donor identities',sex_labels=sorted(obs.sex.astype(str).unique().tolist()),patient_level_inference=False))
table='\n'.join(f'- {r.celltype}: cancer n={r.n},adjacent n={r.n_reference}; means={r.mean_tumor:.6g}/{r.mean_adjacent:.6g}; detection={r.detected_fraction_tumor:.1%}/{r.detected_fraction_adjacent:.1%}; P={r.p_value:.6g}' for r in d[d.annotation_level.eq('minor')].itertuples())
(pub/'README_CN.md').write_text(f'''# SLC6A6 髓系亚群直接细胞比较

本轮问题：作者是否有可追溯的髓系亚群注释；同亚群肿瘤与癌旁细胞的SLC6A6表达如何？用户后续明确选择直接细胞比较，故本轮不进行供者配对或等权。

输入与范围：原PRAD24固定CELLxGENE版本，髓系共{len(my)}细胞，7个minor类别和13个subset。原有字段与作者元数据字典、20221004注释脚本和20230714更新代码一致，包含mature_DC更名mDC_regs、CLEC10A/FCGR3A/CLEC9A/IRF4等更新。当前沿用作者v2注释，不重跑聚类，不以单一标志基因重新定类。来源代码固定提交 {a.author_commit}，来源/哈希见source_manifest.tsv。

实际结果：
{table}

新手解释：每类细胞分别合并全部肿瘤细胞与全部癌旁细胞，展示平均表达和检测比例，Mann-Whitney检验比较分布。P<0.05按细胞层面标记，q值保留。不同患者贡献细胞数不等，P不是患者层面的证据。minor主比较与subset探索比较分别列族，不只挑显著项。

限制/反证：样本仍来自同一研究；癌旁不是健康人。作者标签可追溯，但不是独立金标准验证，特别是少量DC亚群。标志基因图仅作表达一致性检查，不是新注释或功能证据。作者泛化ontology中的部分DC仅标为myeloid cell；本轮使用更细的作者标签而非将ontology强行细化。此前每侧20细胞、至少8供者配对标准下，细分亚群覆盖不足；本轮分析单位依用户要求改为细胞，不声称克服了供者覆盖限制，也不覆盖既往供者分析。表达不等于转运活性。

当前决定：仅用于直接表达对照与假设生成，保留所有结果及细胞数。原有供者配对报告不修改。

下一步：若进行患者层面的验证，仍需更多配对供者/独立数据。

复现：python prad_slc6a6_myeloid_cells_v1.py --out NEW_RUN --source-run ORIGINAL_RUN --code-commit COMMIT --author-commit AUTHOR_COMMIT。源矩阵、细胞值、细胞身份留server165；公开汇总与图。
''',encoding='utf-8')
manifest=[dict(source_id=f.name,path_or_url=str(f),sha256=sha(f)) for f in [src,lp,Path(__file__),a.source_run/'private/sc_donor_profiles_private.tsv']]
for f in sorted((out/'source/author_annotation').iterdir()):
 path=f.name.replace('__','/');manifest.append(dict(source_id=f.name,path_or_url=f'https://raw.githubusercontent.com/swarbricklab/apostolov_pca_atlas/{a.author_commit}/{path}',sha256=sha(f)))
save(pd.DataFrame(manifest),pub/'source_manifest.tsv')
save(pd.DataFrame([dict(file=f.name,sha256=sha(f)) for f in sorted(pub.iterdir())]),pub/'checksums.tsv')
adata.file.close();(out/'.running').rename(out/'COMPLETE')
print(table,flush=True)
print(d[d.annotation_level.eq('subset')][['celltype','n','n_reference','mean_tumor','mean_adjacent','p_value']].to_string(index=False))
