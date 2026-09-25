"""Donor-equal descriptive cell sources; no expression hypothesis tests."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np,pandas as pd
from gbm_discovery_v1 import save,js,sha,check
V='gbm_sc_source_v1b'
MAP={'astrocyte':'Astrocyte','oligodendrocyte':'Oligodendrocyte','neuron':'Neuron','myeloid cell':'Myeloid','neoplastic cell':'Neoplastic','oligodendrocyte precursor cell':'OPC','vascular lymphangioblast':'Vascular','malignant cell':'Neoplastic','macrophage':'Myeloid','microglial cell':'Myeloid','dendritic cell':'Myeloid','monocyte':'Myeloid','mast cell':'Myeloid','mature T cell':'Lymphoid','plasma cell':'Lymphoid','natural killer cell':'Lymphoid','B cell':'Lymphoid','mural cell':'Vascular','endothelial cell':'Vascular','radial glial cell':'Other glial'}
def main(out,commit):
 pub=out/'public/06_EXTERNAL_v1b';pub.mkdir(exist_ok=False);genes=pd.read_csv(out/'source/genes_unique.tsv',sep='\t');profiles=[];ranks=[];coverage=[];maps=[];registry=[];private=[]
 spec=dict(version=V,code_commit=commit,minimum_cells_per_donor_category=20,minimum_distinct_source_donors_per_category=3,minimum_categories_for_ranking=2,minimum_detection_for_rank=.01,bootstrap=1000,master_seed=20260925,bootstrap_unit='whole source donor;duplicate draws are weights,not new donors;at least3 distinct eligible original labels per category in each draw',normalization='log1p(10000*counts/full_retained_raw_gene_library);not UMI for Smartseq',gene_match='exact unique gene symbol to source Ensembl;missing kept',cell_annotation='Darmanis author classes standardized by curator;Neftel GBmap reannotation',scope='all candidate142 genes;coarse category descriptive sources;original tissue regions retained and flagged;no subtype mechanism inference',P_q='NA;no tests',known_overlap_removed='only author=Neftel2019 selected from GBmap;Darmanis2017 part of GBmap excluded',independence_limit='different source studies;cross-publication patient alias equivalence not independently verified')
 js(pub/'analysis_spec.json',spec)
 for study in ['Darmanis2017','Neftel2019']:
  summary=json.loads((out/'source'/f'{study}_extraction_summary.json').read_text());registry.append(summary);npz=np.load(out/'private'/f'{study}_candidate_expression.npz');z=npz['z'];det=npz['detected'];measured=npz['measured'];obs=pd.read_csv(out/'private'/f'{study}_cell_metadata.tsv',sep='\t');assert obs.donor_id.notna().all() and not obs.donor_id.isin(['unknown','']).any();assert set(obs.cell_type)<=set(MAP);obs['coarse']=obs.cell_type.map(MAP);donors=sorted(obs.donor_id.unique());types=sorted(obs.coarse.unique());nd=len(donors);nc=len(types);ng=len(genes)
  means=np.full((nd,nc,ng),np.nan);detections=means.copy();counts=np.zeros((nd,nc),int)
  for i,donor in enumerate(donors):
   for j,ct in enumerate(types):
    use=obs.donor_id.eq(donor)&obs.coarse.eq(ct);n=int(use.sum());counts[i,j]=n
    if n>=20:
     means[i,j,measured]=z[use][:,measured].mean(0);detections[i,j,measured]=det[use][:,measured].mean(0)
     for g,k in enumerate(genes.gene):private.append(dict(study=study,donor_id=donor,celltype=ct,gene=k,n_cells=n,mean_expression=means[i,j,g],detection=detections[i,j,g]))
  valid=np.sum(counts>=20,axis=0)>=3;avg=np.divide(np.nansum(means,axis=0),np.sum(np.isfinite(means),axis=0),out=np.full((nc,ng),np.nan),where=np.sum(np.isfinite(means),axis=0)>0);avg[~valid]=np.nan
  detection=np.divide(np.nansum(detections,axis=0),np.sum(np.isfinite(detections),axis=0),out=np.full((nc,ng),np.nan),where=np.sum(np.isfinite(detections),axis=0)>0);detection[~valid]=np.nan
  for j,ct in enumerate(types):
   for g,k in enumerate(genes.gene):
    vals=means[:,j,g];vals=vals[np.isfinite(vals)];ok=bool(valid[j] and measured[g]);profiles.append(dict(cancer='GBM',cohort=study,stage_id='06_EXTERNAL',run_id=out.name,analysis_version=V,analysis_type='donor_equal_cell_source',metabolite_key='NA',metabolite_name='NA',gene=k,unit='source_donor_label',n=int(len(vals)),n_reference=np.nan,effect_type='donor_equal_mean_log1pCP10k',effect=avg[j,g],ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,q_value=np.nan,test_family='SC_SOURCE_DESCRIPTIVE',family_n_evaluable=np.nan,status='DONE' if ok else 'NOT_EVALUABLE',reason='NA' if ok else 'gene_missing_or_fewer_than3_donors_with20_cells',source_id=study,stable_gene_id=genes.stable_gene_id.iloc[g],partition='all_original_regions',celltype=ct,n_source_labels_total=nd,n_source_labels_eligible=len(vals),n_cells_total=int(counts[:,j].sum()),n_cells_eligible=int(counts[counts[:,j]>=20,j].sum()),mean_detection_fraction=detection[j,g],median_donor_mean=float(np.median(vals)) if ok else np.nan,q25=float(np.quantile(vals,.25)) if ok else np.nan,q75=float(np.quantile(vals,.75)) if ok else np.nan,annotation_origin=summary['annotation_origin']))
  # Common whole-donor draws shared across genes, without creating new independent labels.
  rng=np.random.default_rng(int(hashlib.sha256((V+study).encode()).hexdigest()[:8],16));draws=rng.integers(nd,size=(1000,nd));weights=np.stack([np.bincount(a,minlength=nd) for a in draws]);wins=np.zeros((nc,ng));validboot=np.zeros(ng,int)
  for w in weights:
   eligible=(counts>=20)&(w[:,None]>0);catok=valid&(eligible.sum(0)>=3);den=(np.isfinite(means)*w[:,None,None]).sum(0);av=np.divide((np.nan_to_num(means)*w[:,None,None]).sum(0),den,out=np.full((nc,ng),np.nan),where=den>0);av[~catok]=np.nan
   for g in np.flatnonzero(measured):
    good=np.isfinite(av[:,g])
    if good.sum()<2:continue
    top=np.nanmax(av[:,g]);tie=good&np.isclose(av[:,g],top,atol=1e-10,rtol=0);wins[tie,g]+=1/tie.sum();validboot[g]+=1
  for g,k in enumerate(genes.gene):
   ix=np.flatnonzero(np.isfinite(avg[:,g]));row=dict(study=study,gene=k,stable_gene_id=genes.stable_gene_id.iloc[g],status='NOT_EVALUABLE',reason='fewer_than2_eligible_categories_or_gene_missing',top_celltype='NA',runner_celltype='NA',top_gap=np.nan,bootstrap_top_frequency=np.nan,bootstrap_valid=int(validboot[g]),rank_tie=False)
   if len(ix)>=2 and np.nanmax(detection[:,g])>=.01:
    ordered=ix[np.argsort(-avg[ix,g])];top=avg[ordered[0],g];ties=ix[np.isclose(avg[ix,g],top,rtol=0,atol=1e-10)];row.update(status='DONE',reason='descriptive_not_functional_evidence',top_celltype=';'.join(types[j] for j in ties),runner_celltype=types[ordered[1]],top_gap=float(top-avg[ordered[1],g]),rank_tie=len(ties)>1,bootstrap_top_frequency=float(wins[ties,g].sum()/validboot[g]) if validboot[g] else np.nan)
   elif len(ix)>=2:row['reason']='maximum_mean_detection_below1percent'
   ranks.append(row)
  coverage.append(pd.read_csv(out/'source'/f'{study}_gene_coverage.tsv',sep='\t'));maps += [dict(study=study,original_celltype=ct,coarse_celltype=MAP[ct],annotation_origin=summary['annotation_origin']) for ct in sorted(obs.cell_type.unique())]
  print('summarized',study,flush=True)
 save(pd.DataFrame(private),out/'private/sc_donor_profiles_private.tsv');prof=pd.DataFrame(profiles);rank=pd.DataFrame(ranks);save(prof,pub/'sc_celltype_profiles.tsv');save(rank,pub/'sc_source_stability.tsv');save(pd.concat(coverage),pub/'sc_gene_coverage.tsv');save(pd.DataFrame(registry),pub/'sc_dataset_registry.tsv');save(pd.DataFrame(maps),pub/'sc_annotation_map.tsv')
 cross=[]
 for _,g in genes.iterrows():
  a=rank[(rank.study=='Darmanis2017')&(rank.gene==g.gene)].iloc[0];b=rank[(rank.study=='Neftel2019')&(rank.gene==g.gene)].iloc[0];row=dict(gene=g.gene,stable_gene_id=g.stable_gene_id,study1_top=a.top_celltype,study2_top=b.top_celltype,study1_bootstrap=a.bootstrap_top_frequency,study2_bootstrap=b.bootstrap_top_frequency,same_top_all_categories=np.nan,same_top_shared_categories=np.nan,same_top_both_bootstrap_ge080=np.nan,status='NOT_EVALUABLE')
  if a.status=='DONE' and b.status=='DONE' and not a.rank_tie and not b.rank_tie:row.update(status='DONE',same_top_all_categories=a.top_celltype==b.top_celltype,same_top_both_bootstrap_ge080=bool(a.top_celltype==b.top_celltype and min(a.bootstrap_top_frequency,b.bootstrap_top_frequency)>=.8))
  aa=prof[(prof.cohort=='Darmanis2017')&(prof.gene==g.gene)&(prof.status=='DONE')].set_index('celltype');bb=prof[(prof.cohort=='Neftel2019')&(prof.gene==g.gene)&(prof.status=='DONE')].set_index('celltype');shared=sorted(set(aa.index)&set(bb.index));row['n_shared_categories']=len(shared)
  if len(shared)>=2 and aa.loc[shared,'mean_detection_fraction'].max()>=.01 and bb.loc[shared,'mean_detection_fraction'].max()>=.01:
   av=aa.loc[shared,'effect'];bv=bb.loc[shared,'effect'];at=av.index[np.isclose(av,av.max(),rtol=0,atol=1e-10)].tolist();bt=bv.index[np.isclose(bv,bv.max(),rtol=0,atol=1e-10)].tolist();row['same_top_shared_categories']=at==bt if len(at)==len(bt)==1 else np.nan
  cross.append(row)
 save(pd.DataFrame(cross),pub/'sc_cross_study.tsv');summary=dict(status='DONE',candidate_genes=len(genes),studies=registry,eligible_profile_rows=int(prof.status.eq('DONE').sum()),evaluable_rank_rows=int(rank.status.eq('DONE').sum()),same_top_all=int(pd.DataFrame(cross).same_top_all_categories.eq(True).sum()),shared_top_same=int(pd.DataFrame(cross).same_top_shared_categories.eq(True).sum()),no_P_q_tests=True,full_gene_library_normalization=True,donor_equal_weighting=True,all_genes_retained=True);js(pub/'validation.json',summary)
 paths=[out/'source/genes_unique.tsv',out/'source/Darmanis2017_cellxgene.h5ad',Path(__file__)]+[out/'private'/f'{s}_candidate_expression.npz' for s in ['Darmanis2017','Neftel2019']]+[out/'private'/f'{s}_cell_metadata.tsv' for s in ['Darmanis2017','Neftel2019']];save(pd.DataFrame([dict(source_id=p.name,path_or_url=str(p),sha256=sha(p)) for p in paths]),pub/'source_manifest.tsv')
 (pub/'README_CN.md').write_text('# GBM全部候选的两研究细胞表达来源\n\n## 本轮问题\n142候选基因主要在哪些细胞类别表达？\n\n## 输入与范围\nDarmanis2017原始3589细胞/4供者标签；GBmap内仅Neftel2019子集19030细胞/23供者标签，排除图谱中Darmanis重复数据。沿用作者/图谱注释，不重聚类。\n\n## 实际结果\n'+json.dumps({k:v for k,v in summary.items() if k!='studies'},indent=2)+'\n\n## 新手解释\n先在每供者每类别汇总，再等权平均。至少20细胞、3来源供者才展示；灰色/NA表示不可评估。1000次整供者重采样的最高类别保持率仅为描述。\n\n## 限制/反证\n两研究注释粒度、组织取样与归一化背景不同；全部原取样区域一起描述，不称纯肿瘤核心。跨研究供者别名未独立核验，不称代谢轴独立验证。最高表达不等于唯一作用细胞；没有酶活、通量或机制推断。\n\n## 当前决定\n保留全基因及缺测记录，按共同粗类别和各研究全部类别分别比较。\n\n## 下一步\n回接关系级/基因级比较表，生成全候选点图与热图。\n\n## 复现\n先运行gbm_extract_sc_v1.py，再运行gbm_sc_source_v1.py --out AUDITED_RUN --code-commit COMMIT。\n',encoding='utf-8');check(pub);print(json.dumps({k:v for k,v in summary.items() if k!='studies'}),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--code-commit',required=True);a=p.parse_args();main(a.out,a.code_commit)
