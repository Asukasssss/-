"""GSE159115 whole-library normalization and donor-equal source descriptions."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import argparse,json,hashlib,platform
from pathlib import Path
import numpy as np,pandas as pd,h5py
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ccrcc_stats_v1 import PREFIX,save,js,sha,checksums
R=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716');S=R/'data/candidates/gse159115_kyn_cell_source_v0.1';V='ccrcc_sc_source_v1'
def main(out,commit):
 with (out/'.sc_running').open('x') as f:f.write(V)
 pub=out/'public/06_EXTERNAL';pub.mkdir(exist_ok=False)
 gp=out/'source/gene_membership.tsv';genes=pd.read_csv(gp,sep='\t').sort_values('gene');assert genes.gene.is_unique
 ap=S/'GSE159115_ccRCC_anno.csv.gz';obs=pd.read_csv(ap);assert obs.cell.is_unique and obs.patient.notna().all();obs=obs.set_index('cell')
 spec=dict(version=V,code_commit=commit,study='GSE159115',scope='author ccRCC tumor annotation only;normal excluded from primary source;all current/historical genes',
  annotation='author anno unchanged;ua retained as unassigned and not ranked;no reclustering or inferred malignancy',
  donor='author patient field,all samples of same label aggregated',normalization='log1p(10000*integer_counts/sum_all_H5_gene_counts_per_cell)',
  detection='counts>0',gene_match='exact stable Ensembl when available,otherwise exact unique primary symbol;no inferred aliases',
  minimum_cells_per_donor_type=20,minimum_donors_per_type=3,ranking='>=2 evaluable assigned categories;max detection>=0.01;ties split credit',
  bootstrap='1000 whole donor draws, all types together;each bootstrap category requires>=3 distinct eligible original donors',seed=20260925,
  plot='alphabetical genes and author cell labels;no clustering;row z-score heatmap display only, symmetric[-2,2];raw donor-equal mean for dot color;missing grey',
  independent_studies=1,P_q='NA descriptive',software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,h5py=h5py.__version__,matplotlib=matplotlib.__version__))
 js(pub/'analysis_spec.json',spec)
 norm=pd.DataFrame(np.nan,index=obs.index,columns=genes.gene);det=norm.copy();lib=pd.Series(np.nan,index=obs.index);seen=set();coverage=[];used=[];gene_counts=set()
 for path in sorted((S/'h5').glob('*.h5')):
  sample='SI_'+path.name.split('_SI_',1)[1].split('_filtered')[0] if '_SI_' in path.name else path.name.split('_filtered')[0]
  if sample not in set(obs['sample']):continue
  with h5py.File(path,'r') as f:
   z=f['GRCh38'];symbols=np.array([v.decode() for v in z['gene_names'][:]]);ids=np.array([v.decode().split('.')[0] for v in z['genes'][:]]);bars=[sample+'_'+v.decode() for v in z['barcodes'][:]]
   x=sparse.csc_matrix((z['data'][:],z['indices'][:],z['indptr'][:]),shape=(len(symbols),len(bars)))
   assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all();gene_counts.add(len(symbols))
   j=np.flatnonzero(pd.Index(bars).isin(obs.index));cells=np.array(bars)[j];assert not set(cells)&seen;seen.update(cells);x=x[:,j]
   library=np.asarray(x.sum(0)).ravel();assert (library>0).all();lib.loc[cells]=library
   for _,g in genes.iterrows():
    hit=np.flatnonzero(ids==g.stable_gene_id) if str(g.stable_gene_id).startswith('ENSG') else np.flatnonzero(symbols==g.gene)
    coverage.append(dict(source_file=path.name,gene=g.gene,stable_gene_id=g.stable_gene_id,n_hits=len(hit),status='DONE' if len(hit)==1 else 'NOT_EVALUABLE',reason='exact_unique_match' if len(hit)==1 else 'missing_or_ambiguous_feature'))
    if len(hit)==1:
     values=x[hit[0],:].toarray().ravel();norm.loc[cells,g.gene]=np.log1p(10000*values/library);det.loc[cells,g.gene]=(values>0).astype(float)
   used.append(path)
  print('SC_FILE',len(used),'matched_cells',len(seen),flush=True)
 assert seen==set(obs.index),f'Unmatched annotated cells {len(set(obs.index)-seen)}'
 save(pd.DataFrame(coverage),pub/'sc_feature_coverage.tsv');save(pd.DataFrame({'library':lib}),out/'private/sc_library_private.tsv')
 profiles=[]
 for (donor,ct),z in obs.groupby(['patient','anno']):
  for _,g in genes.iterrows():
   vals=norm.loc[z.index,g.gene];dv=det.loc[z.index,g.gene]
   profiles.append(dict(donor=donor,celltype=ct,gene=g.gene,n_cells=len(z),n_cells_measured=int(vals.notna().sum()),mean_expression=float(vals.mean()) if vals.notna().all() else np.nan,detection_fraction=float(dv.mean()) if dv.notna().all() else np.nan))
 prof=pd.DataFrame(profiles);save(prof,out/'private/sc_donor_profiles_private.tsv')
 donors=sorted(obs.patient.unique());cats=sorted(obs.anno.unique());rng=np.random.default_rng(20260925);draws=rng.integers(len(donors),size=(1000,len(donors)));weights=np.array([np.bincount(v,minlength=len(donors)) for v in draws]);public=[];ranks=[]
 for _,g in genes.iterrows():
  pg=prof[prof.gene.eq(g.gene)];eligible=pg[pg.n_cells.ge(20)&pg.mean_expression.notna()];values={};detect={}
  for ct in cats:
   z=eligible[eligible.celltype.eq(ct)];allz=pg[pg.celltype.eq(ct)];okay=len(z)>=3 and ct!='ua'
   d=dict.fromkeys(PREFIX,np.nan);d.update(cancer='ccRCC',cohort='GSE159115',stage_id='06_EXTERNAL',run_id=out.name,analysis_version=V,analysis_type='sc_expression_source',gene=g.gene,stable_gene_id=g.stable_gene_id,metabolite_key='NA',metabolite_name='NA',unit='author_source_donor_label',n=len(z),n_reference=len(donors),effect_type='donor_equal_mean_log1p_CP10K',test_family='SC_SOURCE_DESCRIPTIVE',status='DONE' if okay else 'NOT_EVALUABLE',reason='descriptive_only' if okay else 'unassigned_label' if ct=='ua' else 'insufficient_coverage_or_gene_missing',source_id='GSE159115',celltype=ct,partition='tumor',n_cells_total=int(allz.n_cells.sum()),n_cells_eligible=int(z.n_cells.sum()),annotation_origin='author_anno',normalization=spec['normalization'])
   if okay:
    d.update(effect=z.mean_expression.mean(),mean_detection_fraction=z.detection_fraction.mean(),median_donor_mean=z.mean_expression.median(),q25=z.mean_expression.quantile(.25),q75=z.mean_expression.quantile(.75));values[ct]=d['effect'];detect[ct]=d['mean_detection_fraction']
   public.append(d)
  rr=dict(gene=g.gene,stable_gene_id=g.stable_gene_id,top_celltype='NA',runner_celltype='NA',top_gap=np.nan,bootstrap_top_frequency=np.nan,bootstrap_valid=0,n_eligible_celltypes=len(values),status='NOT_EVALUABLE',reason='fewer_than2_types_or_max_detection_below1pct',study='GSE159115')
  if len(values)>=2 and max(detect.values())>=.01:
   order=sorted(values,key=lambda k:(-values[k],k));tops=[k for k in order if np.isclose(values[k],values[order[0]],rtol=1e-12,atol=1e-12)];cols=list(values)
   mm=eligible.pivot(index='donor',columns='celltype',values='mean_expression').reindex(index=donors,columns=cols).to_numpy();dd=eligible.pivot(index='donor',columns='celltype',values='detection_fraction').reindex(index=donors,columns=cols).to_numpy();valid=np.isfinite(mm);den=weights@valid;unique=(weights>0).astype(int)@valid.astype(int)
   with np.errstate(invalid='ignore',divide='ignore'):bm=(weights@np.nan_to_num(mm))/den;bd=(weights@np.nan_to_num(dd))/den
   bm[unique<3]=np.nan;bd[unique<3]=np.nan;ok=(np.isfinite(bm).sum(1)>=2)&(np.max(np.where(np.isfinite(bd),bd,-1),axis=1)>=.01)
   ties=np.isclose(bm[ok],np.nanmax(bm[ok],axis=1)[:,None],rtol=1e-12,atol=1e-12);credit=ties/ties.sum(1)[:,None]
   freq=float(credit[:,[cols.index(c) for c in tops]].sum(1).mean()) if ok.any() else np.nan
   rr.update(top_celltype=';'.join(tops),runner_celltype=order[1],top_gap=values[order[0]]-values[order[1]],bootstrap_top_frequency=freq,bootstrap_valid=int(ok.sum()),status='DONE',reason='top_tied_descriptive' if len(tops)>1 else 'descriptive_only')
  ranks.append(rr)
 d=pd.DataFrame(public);d=d[PREFIX+[c for c in d if c not in PREFIX]];rank=pd.DataFrame(ranks);save(d,pub/'sc_celltype_profiles.tsv');save(rank,pub/'sc_source_stability.tsv')
 save(pd.DataFrame({'author_label':cats,'display_label':cats,'annotation_origin':'author_anno_unchanged'}),pub/'sc_annotation_map.tsv')
 matrix=d.pivot(index='gene',columns='celltype',values='effect').reindex(genes.gene);values=matrix.to_numpy();mean=np.nanmean(values,1,keepdims=True);sd=np.nanstd(values,1,keepdims=True);z=(values-mean)/np.where(sd>0,sd,1)
 cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#d9d9d9')
 fig,ax=plt.subplots(figsize=(12,max(12,.22*len(genes))));im=ax.imshow(z,aspect='auto',interpolation='nearest',cmap=cmap,vmin=-2,vmax=2);ax.set_yticks(range(len(genes)));ax.set_yticklabels(matrix.index,fontsize=7);ax.set_xticks(range(len(matrix.columns)));ax.set_xticklabels(matrix.columns,rotation=45,ha='right');ax.set_title('ccRCC all-candidate cell expression source | GSE159115\nDonor-equal; row z-score for display only; grey = not evaluable');fig.colorbar(im,ax=ax,label='Row z-score');fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(pub/('all_gene_source_heatmap.'+ext),dpi=180)
 plt.close(fig)
 detection=d.pivot(index='gene',columns='celltype',values='mean_detection_fraction').reindex(index=matrix.index,columns=matrix.columns).to_numpy();xx,yy=np.meshgrid(np.arange(len(matrix.columns)),np.arange(len(matrix)));valid=np.isfinite(values)
 fig,ax=plt.subplots(figsize=(12,max(12,.22*len(genes))));ax.scatter(xx[~valid],yy[~valid],s=12,c='#d9d9d9',marker='s');im=ax.scatter(xx[valid],yy[valid],s=90*detection[valid],c=values[valid],cmap='viridis',vmin=0);ax.set_yticks(range(len(genes)));ax.set_yticklabels(matrix.index,fontsize=7);ax.set_xticks(range(len(matrix.columns)));ax.set_xticklabels(matrix.columns,rotation=45,ha='right');ax.invert_yaxis();ax.set_title('ccRCC donor-equal expression | GSE159115\nDot area = detection fraction; colour = log1p(CP10K) mean');fig.colorbar(im,ax=ax,label='Donor-equal mean');fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(pub/('all_gene_source_dotplot.'+ext),dpi=180)
 plt.close(fig)
 cov=pd.DataFrame(coverage);val=dict(status='DONE',cells=len(obs),author_donor_labels=len(donors),genes=len(genes),genes_measured_all_files=int(cov.groupby('gene').status.apply(lambda z:z.eq('DONE').all()).sum()),source_rank_evaluable=int(rank.status.eq('DONE').sum()),all_annotated_cells_matched=True,raw_counts_nonnegative_integer=True,full_gene_denominators=sorted(gene_counts),donor_equal=True,second_independent_study_completed=False,new_UMAP=False,new_clustering=False)
 js(pub/'validation.json',val)
 save(pd.DataFrame([dict(source_id='GSE159115',status='DONE',source_url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159115',scope='author ccRCC tumor cells',limitation='single study;source labels not genotype reverified'),dict(source_id='GSE269826',status='NOT_RUN',source_url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269826',scope='second study candidate',limitation='existing downloads have aria2 partial markers;not analyzed or counted as replication')]),pub/'sc_dataset_registry.tsv')
 save(pd.DataFrame([dict(source_id=p.name,path_or_url=str(p),sha256=sha(p),access_scope='SERVER_PRIVATE' if p in used or p==ap else 'CODE_OR_SUMMARY') for p in [ap,gp,Path(__file__),*used]]),pub/'source_manifest.tsv')
 (pub/'README_CN.md').write_text(f'''# ccRCC全候选单细胞表达来源\n\n问题：当前/历史并集基因主要在哪些作者细胞类别表达？\n\n输入：GSE159115作者ccRCC肿瘤注释及10x计数，{len(obs)}细胞、{len(donors)}个作者患者标签、{len(genes)}个候选基因。\n\n实际结果：{val['genes_measured_all_files']}基因所有文件精确覆盖，{val['source_rank_evaluable']}基因可描述最高表达类别；完整点图和热图随附。\n\n新手解释：先供者内再供者间等权，颜色并不是酶活；最高表达类别不能确定代谢物来源。\n\n限制：目前一项研究，不能称双研究复现。低覆盖为NA；ua不参与排名。第二资料下载不完整。作者注释没有UMAP坐标，本轮不新算UMAP。\n\n当前决定：仅作为表达背景，和患者统计分列整合。\n\n下一步：关系/基因表与报告。\n\n复现：python ccrcc_sc_source_v1.py --out RUN --code-commit COMMIT；每供者类别至少20细胞、每类至少3供者，1000次联合供者bootstrap。\n''',encoding='utf-8')
 checksums(pub);(out/'.sc_running').rename(out/'.sc_done');print(json.dumps(val),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--code-commit',required=True);a=p.parse_args();main(a.out,a.code_commit)
