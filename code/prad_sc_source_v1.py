"""Server-only donor-equal source profiles; preserve author labels, no reclustering."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import argparse,json,hashlib,platform
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
V='prad_sc_source_v1'
PREFIX='cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for z in iter(lambda:f.read(1048576),b''):h.update(z)
 return h.hexdigest()
def js(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
def main(out,commit):
 assert (out/'.running').read_text()=='prad_discovery_v1'
 pub=out/'public/06_EXTERNAL';pub.mkdir(exist_ok=False)
 src=out/'source/PRAD24_cellxgene.h5ad';edgesp=out/'source/direct_relations.tsv';genes=pd.read_csv(edgesp,sep='\t')[['gene','stable_gene_id']].drop_duplicates().sort_values('gene');assert genes.gene.is_unique
 a=ad.read_h5ad(src,backed='r');obs=a.obs.copy();assert obs.index.is_unique and obs.donor_id.notna().all()
 assert a.raw is not None and a.raw.var.index.is_unique
 assert obs.groupby('sample',observed=True).donor_id.nunique().max()==1
 obs['source_celltype']=obs.celltype_major_v2.astype(str)
 epithelial=obs.source_celltype.eq('Epithelial')
 obs.loc[epithelial,'source_celltype']='Epithelial_'+obs.loc[epithelial,'malignant_anno_merged'].astype(str)
 assert obs['type'].isin(['cancer','adj_benign']).all()
 ann=obs[['celltype_major_v2','malignant_anno_merged','source_celltype']].drop_duplicates()
 ann['annotation_origin']='Depositor author-v2 annotation via CELLxGENE; no inference/reclustering';save(ann,pub/'sc_annotation_map.tsv')
 spec=dict(version=V,code_commit=commit,dataset='PRAD24_Apostolov_CELLxGENE',dataset_version='68b23fda-7191-46a5-8870-819feca3e66e',
  source_publication='10.1101/2024.10.23.619925',publication_status='preprint per EuropePMC record on2026-09-22',
  accession_warning='Collection links GSE145843, which describes normal prostate elsewhere; do not adopt that accession. Use exact CELLxGENE dataset version and matching embedded title/citation; accession reconciliation pending',
  primary_partition='type=cancer;not all cells carry malignant annotation',supportive_partition='type=adj_benign;same study and overlapping donors,not independent replication',
  donor_field='donor_id;merge multiple sample blocks within donor and tissue partition',label='celltype_major_v2;Epithelial split by depositor malignant_anno_merged',
  count_location='raw.X',normalization='log1p(10000*raw_count/sum_all_raw_genes_per_cell)',detection='raw_count>0',
  per_donor_min_cells=20,per_celltype_min_distinct_donors=3,donor_weight='equal',
  candidate_genes=len(genes),gene_match='exact unique raw.var.feature_name;HGNC stable dictionary from biochemical stage;no synonym guesses',
  bootstrap=1000,seed=20260922,bootstrap_unit='whole donor;all cell types together;distinct original eligible labels>=3 for each bootstrap category',
  ranking='at least2 eligible classes,max mean detection>=1%;Unassigned excluded;ties split credit;80pct descriptive only',
  scientific_scope='expression source only;no cell communication,trajectory,networks,mechanism,subtype significance,or new UMAP',
  P_q='NA;SC_SOURCE_DESCRIPTIVE',software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,anndata=ad.__version__))
 js(pub/'analysis_spec.json',spec)
 var=a.raw.var;symbol=var.feature_name.astype(str);coverage=[];ix=[];measured=[]
 for _,g in genes.iterrows():
  hits=np.flatnonzero(symbol.to_numpy()==g.gene)
  status='DONE' if len(hits)==1 else 'NOT_EVALUABLE'
  coverage.append(dict(gene=g.gene,stable_gene_id=g.stable_gene_id,status=status,reason='exact_unique_gene_symbol' if len(hits)==1 else 'gene_not_present' if len(hits)==0 else 'ambiguous_multiple_features',n_feature_hits=len(hits),ensembl_id=str(var.index[hits[0]]) if len(hits)==1 else 'NA'))
  if len(hits)==1:ix.append(hits[0]);measured.append(g.gene)
 save(pd.DataFrame(coverage),pub/'sc_gene_coverage.tsv')
 norm=np.empty((len(obs),len(ix)),dtype=np.float32);det=np.empty((len(obs),len(ix)),dtype=bool);library=np.empty(len(obs));integer_maxerr=0.
 for start in range(0,len(obs),2000):
  end=min(start+2000,len(obs));x=a.raw.X[start:end].tocsr();assert np.isfinite(x.data).all() and (x.data>=0).all()
  integer_maxerr=max(integer_maxerr,float(np.max(abs(x.data-np.round(x.data)))) if x.nnz else 0.)
  assert integer_maxerr==0, 'raw.X not integer counts; do not re-normalize unknown input'
  lib=np.asarray(x.sum(1)).ravel();assert (lib>0).all();library[start:end]=lib
  v=x[:,ix].toarray();det[start:end]=v>0;norm[start:end]=np.log1p(10000*v/lib[:,None])
 print('counts_validated',len(obs),'genes_measured',len(measured),flush=True)
 pd.DataFrame({'library':library},index=obs.index).to_csv(out/'private/sc_library_private.tsv',sep='\t')
 private=[]
 for (partition,donor,ct),idx in obs.groupby(['type','donor_id','source_celltype'],observed=True).indices.items():
  means=norm[idx].mean(0);detection=det[idx].mean(0)
  for j,g in enumerate(measured):private.append(dict(partition=partition,donor=donor,celltype=ct,gene=g,n_cells=len(idx),mean_log1p_cp10k=float(means[j]),detection_fraction=float(detection[j])))
 prof=pd.DataFrame(private);save(prof,out/'private/sc_donor_profiles_private.tsv')
 public=[];rank=[]
 geneid=genes.set_index('gene').stable_gene_id.to_dict()
 for partition,pp in prof.groupby('partition'):
  donor_names=sorted(pp.donor.unique());categories=sorted(pp.celltype.unique())
  rs=np.random.default_rng(int(hashlib.sha256((V+'|'+partition).encode()).hexdigest()[:8],16))
  draw=rs.integers(len(donor_names),size=(1000,len(donor_names)))
  weights=np.array([np.bincount(z,minlength=len(donor_names)) for z in draw])
  for g in genes.gene:
   pg=pp[pp.gene.eq(g)];eligible=pg[pg.n_cells.ge(20)];rowvalues={};rowdetection={}
   for ct in categories:
    z=eligible[eligible.celltype.eq(ct)];allz=pg[pg.celltype.eq(ct)]
    evaluable=len(z)>=3 and ct!='Unassigned'
    d=dict.fromkeys(PREFIX,np.nan);d.update(cancer='PRAD',cohort='PRAD24_CELLxGENE',stage_id='06_EXTERNAL',run_id=out.name,analysis_version=V,analysis_type='donor_equal_expression_source',metabolite_key='NA',metabolite_name='NA',gene=g,unit='source_donor_label',n=len(z),n_reference=len(donor_names),effect_type='donor_equal_mean_log1p_CP10K',test_family='SC_SOURCE_DESCRIPTIVE',status='DONE' if evaluable else 'NOT_EVALUABLE',reason='NA' if evaluable else 'unassigned_annotation' if ct=='Unassigned' else 'gene_not_measured' if not len(pg) else 'fewer_than3_donors_with20cells',source_id=spec['dataset'],stable_gene_id=geneid[g],partition=partition,celltype=ct,
      n_source_labels_total=len(donor_names),n_source_labels_eligible=len(z),n_cells_total=int(allz.n_cells.sum()),n_cells_eligible=int(z.n_cells.sum()),annotation_origin='depositor_author_v2',normalization=spec['normalization'])
    if evaluable:
     d.update(effect=float(z.mean_log1p_cp10k.mean()),mean_detection_fraction=float(z.detection_fraction.mean()),median_donor_mean=float(z.mean_log1p_cp10k.median()),q25=float(z.mean_log1p_cp10k.quantile(.25)),q75=float(z.mean_log1p_cp10k.quantile(.75)))
     rowvalues[ct]=d['effect'];rowdetection[ct]=d['mean_detection_fraction']
    public.append(d)
   rr=dict(gene=g,stable_gene_id=geneid[g],partition=partition,top_celltype='NA',runner_celltype='NA',top_gap=np.nan,bootstrap_top_frequency=np.nan,bootstrap_valid=0,n_eligible_celltypes=len(rowvalues),status='NOT_EVALUABLE',reason='fewer_than2_eligible_celltypes' if len(rowvalues)<2 else 'max_detection_below1pct')
   if len(rowvalues)>=2 and max(rowdetection.values())>=.01:
    order=sorted(rowvalues,key=lambda c:(-rowvalues[c],c));tops=[c for c in order if np.isclose(rowvalues[c],rowvalues[order[0]],atol=1e-12,rtol=1e-12)]
    cols=list(rowvalues)
    mat=eligible.pivot(index='donor',columns='celltype',values='mean_log1p_cp10k').reindex(index=donor_names,columns=cols).to_numpy(float)
    dm=eligible.pivot(index='donor',columns='celltype',values='detection_fraction').reindex(index=donor_names,columns=cols).to_numpy(float)
    valid=np.isfinite(mat);den=weights@valid;unique=(weights>0).astype(int)@valid.astype(int)
    with np.errstate(invalid='ignore',divide='ignore'):
     bmean=(weights@np.nan_to_num(mat))/den;bd=(weights@np.nan_to_num(dm))/den
    bmean[unique<3]=np.nan;bd[unique<3]=np.nan
    ok=(np.isfinite(bmean).sum(1)>=2)&(np.nanmax(bd,axis=1)>=.01)
    bm=bmean[ok];topmax=np.nanmax(bm,axis=1);ties=np.isclose(bm,topmax[:,None],atol=1e-12,rtol=1e-12)
    credit=ties/ties.sum(1)[:,None];freq=float(credit[:,[cols.index(c) for c in tops]].sum(1).mean()) if ok.any() else np.nan
    rr.update(top_celltype=';'.join(tops),runner_celltype=order[1],top_gap=rowvalues[order[0]]-rowvalues[order[1]],bootstrap_top_frequency=freq,bootstrap_valid=int(ok.sum()),status='DONE',reason='top_tied' if len(tops)>1 else 'descriptive_only')
   rank.append(rr)
 d=pd.DataFrame(public);d=d[PREFIX+[c for c in d if c not in PREFIX]];ranks=pd.DataFrame(rank)
 save(d,pub/'sc_celltype_profiles.tsv');save(ranks,pub/'sc_source_stability.tsv')
 # Same-study cancer vs adjacent benign description, never independent replication.
 cross=ranks[ranks.partition.eq('cancer')].merge(ranks[ranks.partition.eq('adj_benign')],on=['gene','stable_gene_id'],suffixes=('_cancer','_adj_benign'),validate='one_to_one')
 cross['same_study_not_independent']=True;save(cross,pub/'sc_within_study_tissue_context.tsv')
 c=d[d.partition.eq('cancer')].pivot(index='gene',columns='celltype',values='effect').reindex(genes.gene)
 c=c.dropna(axis=1,how='all');rawvalues=c.to_numpy(float);mu=np.nanmean(rawvalues,axis=1,keepdims=True);sd=np.nanstd(rawvalues,axis=1,keepdims=True)
 z=(rawvalues-mu)/np.where(sd>0,sd,1)
 cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#d9d9d9')
 fig,ax=plt.subplots(figsize=(12,24));im=ax.imshow(z,aspect='auto',cmap=cmap,vmin=-2,vmax=2,interpolation='nearest');ax.set_yticks(np.arange(len(c)));ax.set_yticklabels(c.index,fontsize=7);ax.set_xticks(np.arange(len(c.columns)));ax.set_xticklabels(c.columns,rotation=45,ha='right',fontsize=9);ax.set_title('PRAD: all candidate genes, donor-equal expression source\nCancer tissue only; row z-score for display; grey = not evaluable',fontsize=12);fig.colorbar(im,ax=ax,label='Row z-score of donor-equal mean');fig.tight_layout();fig.savefig(pub/'all101_source_heatmap.png',dpi=180);fig.savefig(pub/'all101_source_heatmap.pdf');plt.close(fig)
 detection=d[d.partition.eq('cancer')].pivot(index='gene',columns='celltype',values='mean_detection_fraction').reindex(index=c.index,columns=c.columns)
 fig,ax=plt.subplots(figsize=(12,24));xx,yy=np.meshgrid(np.arange(len(c.columns)),np.arange(len(c)));valid=np.isfinite(rawvalues);sc=ax.scatter(xx[valid],yy[valid],s=detection.to_numpy()[valid]*90,c=rawvalues[valid],cmap='viridis',vmin=0);ax.set_yticks(np.arange(len(c)));ax.set_yticklabels(c.index,fontsize=7);ax.set_xticks(np.arange(len(c.columns)));ax.set_xticklabels(c.columns,rotation=45,ha='right',fontsize=9);ax.invert_yaxis();ax.set_title('PRAD cancer tissue: donor-equal expression and detection\nDot area = detection fraction; colour = mean log1p(CP10K)',fontsize=12);fig.colorbar(sc,ax=ax,label='Donor-equal mean log1p(CP10K)');fig.tight_layout();fig.savefig(pub/'all101_source_dotplot.png',dpi=180);fig.savefig(pub/'all101_source_dotplot.pdf');plt.close(fig)
 registry=[dict(source_id=spec['dataset'],status='DONE',role='single_study_cell_expression_context',cells=len(obs),cancer_cells=int(obs['type'].eq('cancer').sum()),cancer_donors=int(obs.loc[obs['type'].eq('cancer'),'donor_id'].nunique()),source_url='https://cellxgene.cziscience.com/collections/bdac7a53-fe34-4f04-8c46-f9bd5297c099',limitation=spec['accession_warning']),
  dict(source_id='GSE141445',status='NEEDS_REVIEW',role='second_tumor_study_candidate',source_url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141445',limitation='Processed expression archive available;author cell annotation and donor identity alignment not yet established;not analyzed'),
  dict(source_id='TISCH2_PRAD',status='ACCESS_BLOCKED',role='possible_second_study_annotation',source_url='https://tisch.compbio.cn/gallery/?cancer=PRAD&species=Human',limitation='Both local proxy and server HTTP endpoints timed out;no downloaded expression/annotation')]
 save(pd.DataFrame(registry),pub/'sc_dataset_registry.tsv')
 val=dict(status='DONE',input_cells=len(obs),primary_cancer_cells=int(obs['type'].eq('cancer').sum()),cancer_donors=24,all_candidate_genes=len(genes),measured_genes=len(measured),missing_genes=[r['gene'] for r in coverage if r['status']!='DONE'],primary_source_rank_evaluable=int(ranks.partition.eq('cancer').mul(ranks.status.eq('DONE')).sum()),raw_integer_count_error=integer_maxerr,full_gene_denominator=len(var),donor_equal=True,new_clustering=False,new_UMAP=False,second_independent_tumor_study_completed=False)
 js(pub/'validation.json',val)
 save(pd.DataFrame([dict(source_id=p.name,path_or_url=str(p),sha256=sha(p)) for p in [src,edgesp,out/'source/PRAD24_download.json',Path(__file__)]]),pub/'source_manifest.tsv')
 (pub/'README_CN.md').write_text(f'''# PRAD 全101候选基因单细胞表达来源\n\n问题：候选基因在肿瘤组织中的哪些细胞类别表达？\n\n输入：Apostolov等24例未治疗局限性PRAD研究的CELLxGENE固定版本。共{len(obs)}细胞，但主分析只用作者type=cancer的{val['primary_cancer_cells']}细胞；adj_benign另列同研究参考，不混入主来源。\n\n实际结果：101基因中{len(measured)}基因精确唯一匹配；缺测{', '.join(val['missing_genes'])}。主肿瘤来源排名可评估{val['primary_source_rank_evaluable']}基因。完整热图、点图和各供者等权汇总随附。\n\n新手解释：先在每个供者内部计算，再给供者等权；不是大样本贡献更多权重。最高表达类别是描述，80%排名保持率不是功能证明，也不是80%患者都有该作用。\n\n限制：目前仅一项肿瘤研究，第二研究注释未就绪，不能称双研究一致。数据集网页GSE145843连接与24例肿瘤标题不一致，本轮只引用精确CELLxGENE版本及嵌入的Apostolov DOI，不使用该可疑GEO号。原记录为预印本；作者v2注释被原样继承，上皮按作者malignant_anno_merged分开，非本轮推断。全部癌旁数据与肿瘤数据有供者重叠，不当独立队列。T-cells作者大类含NK，SMCs大类含周细胞；保持原分组而不宣称纯细胞亚型。\n\n当前决定：仅作细胞表达来源，不推断代谢物生成/消耗细胞或因果；未测基因和低覆盖类别保留NA。\n\n下一步：整合患者结果；第二研究需可靠作者细胞注释及供者映射后再接入。\n\n复现：python prad_sc_source_v1.py --out SAME_RUN --code-commit COMMIT。raw.X全基因库归一化，每供者类别至少20细胞，每类至少3供者；不重跑聚类或UMAP。\n''',encoding='utf-8')
 save(pd.DataFrame([dict(file=p.name,sha256=sha(p)) for p in sorted(pub.iterdir()) if p.name!='checksums.tsv']),pub/'checksums.tsv');print(json.dumps(val),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--code-commit',required=True);a=ap.parse_args();main(a.out,a.code_commit)
