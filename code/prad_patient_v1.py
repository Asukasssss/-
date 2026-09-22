"""Server-only all-direct-relation PRAD associations and paired RNA background."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import argparse,json,hashlib,platform
from pathlib import Path
import numpy as np,pandas as pd,scipy
from scipy import stats
from statsmodels.stats.multitest import multipletests
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
SRC=ROOT/'data/candidates/camp_primary_tissue_multicancer'
V='prad_patient_v1'
PREFIX='cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def js(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def seed(f,k):return int(hashlib.sha256((V+'|20260922|'+f+'|'+k).encode()).hexdigest()[:8],16)
def base(out,f,g):
 d=dict.fromkeys(PREFIX,np.nan);d.update(cancer='PRAD',cohort='CAMP_PRAD',stage_id='03_PATIENT',run_id=out.name,analysis_version=V,analysis_type=f,test_family=f,gene=g,unit='author_case',status='NOT_EVALUABLE',reason='NA',source_id='CAMP_v0.3.4_PRAD',same_cohort_postselection=True);return d
def bh(d):
 for f,ix in d.groupby('test_family').groups.items():
  ok=d.index.isin(ix)&d.p_value.notna();p=d.loc[ok,'p_value'].to_numpy(float);d.loc[ix,'family_n_planned']=len(ix);d.loc[ix,'family_n_evaluable']=len(p)
  if len(p):
   q=multipletests(p,method='fdr_bh')[1];order=np.argsort(p);check=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);assert np.allclose(q[order],check);d.loc[ok,'q_value']=q
 d['evidence_label']=np.select([d.p_value.isna(),d.q_value.lt(.05),d.p_value.lt(.05)],['NOT_EVALUABLE','FDR_SUPPORTED','NOMINAL_EXPLORATORY'],default='NOT_SUPPORTED_THIS_TEST')
 return d[PREFIX+[c for c in d if c not in PREFIX]]
def centered(z,batches=None):
 z=z.copy()
 if batches is None:return z-z.mean(axis=-1,keepdims=True)
 for b in np.unique(batches):
  ix=batches==b;z[...,ix]-=z[...,ix].mean(axis=-1,keepdims=True)
 return z
def assoc(x,y,batches,f,k,partial=False):
 n=len(x);sd=seed(f,k);d=dict(n=n,n_reference=np.nan,n_complete=n,n_unique_author_cases=n,seed=sd,B=9999,p_method='9999_global_case_permutations_plus1',ci_method='4000_whole_case_bootstrap_percentile',ci_estimand='Spearman_rho',effect_type='Spearman_rho',bootstrap_valid=0)
 if n<8:return dict(d,status='NOT_EVALUABLE',reason='fewer_than8_complete_tumors')
 if len(np.unique(x))<2 or len(np.unique(y))<2:return dict(d,status='NOT_EVALUABLE',reason='constant_input')
 rng=np.random.default_rng(sd);rx=stats.rankdata(x);ry=stats.rankdata(y)
 a=centered(rx,batches if partial else None);b=centered(ry,batches if partial else None)
 denom=np.sqrt((a*a).sum()*(b*b).sum())
 if denom<=0:return dict(d,status='NOT_EVALUABLE',reason='zero_residual_rank_variance')
 rho=float(a@b/denom)
 if partial:
  perms=np.tile(np.arange(n),(9999,1))
  for label in np.unique(batches):
   ix=np.flatnonzero(batches==label);perms[:,ix]=ix[np.argsort(rng.random((9999,len(ix))),axis=1)]
  d.update(p_method='9999_within_author_batch_case_permutations_plus1',effect_type='partial_Spearman_rank_residual_r',ci_estimand='partial_Spearman_rank_residual_r',ci_method='4000_batch_stratified_whole_case_bootstrap_percentile')
 else:perms=np.argsort(rng.random((9999,n)),axis=1)
 rp=(b[perms]@a)/denom;p=float((1+(abs(rp)>=abs(rho)-1e-12).sum())/10000)
 if not partial:assert abs(rho-stats.spearmanr(x,y).statistic)<1e-12
 if partial:
  ix=np.tile(np.arange(n),(4000,1))
  for label in np.unique(batches):
   sel=np.flatnonzero(batches==label);ix[:,sel]=sel[rng.integers(len(sel),size=(4000,len(sel)))]
 else:ix=rng.integers(n,size=(4000,n))
 bx=centered(stats.rankdata(x[ix],axis=1),batches if partial else None);by=centered(stats.rankdata(y[ix],axis=1),batches if partial else None)
 den=np.sqrt((bx*bx).sum(1)*(by*by).sum(1));valid=den>0;boot=(bx[valid]*by[valid]).sum(1)/den[valid]
 d.update(effect=rho,rho=rho,p_value=p,status='DONE',reason='NA',bootstrap_valid=int(valid.sum()),n_author_batches=len(np.unique(batches)))
 if valid.sum()>=3600:d.update(ci_lower=float(np.quantile(boot,.025)),ci_upper=float(np.quantile(boot,.975)))
 else:d.update(ci_lower=np.nan,ci_upper=np.nan,reason='bootstrap_valid_fraction_below90pct')
 return d
def main(out,commit):
 assert (out/'.running').read_text()=='prad_discovery_v1'
 pub=out/'public/03_PATIENT';pub.mkdir(exist_ok=False)
 edgesp=out/'source/direct_relations.tsv';edges=pd.read_csv(edgesp,sep='\t');assert not edges.duplicated(['feature_id','stable_gene_id']).any()
 genes=edges[['gene','stable_gene_id']].drop_duplicates();assert genes.gene.is_unique
 mp=out/'private/tumor_multiomics_map_private.tsv';m=pd.read_csv(mp,sep='\t',dtype=str);assert m['case'].is_unique
 pp=out/'private/RNA_pairs_private.tsv';pairs=pd.read_csv(pp,sep='\t',dtype=str)
 rp=SRC/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/m.RNAFile.iloc[0]
 xp=SRC/'processed_metabolomics/PreprocessedData_PRAD.xlsx'
 rna=pd.read_csv(rp,index_col=0);matrix=pd.read_excel(xp,sheet_name='data_imputed',index_col=0);raw=pd.read_excel(xp,sheet_name='data',index_col=0)
 for d in [matrix,raw]:d.index=d.index.astype(str).str.strip();assert d.index.is_unique
 assert rna.index.is_unique and rna.columns.is_unique
 spec=dict(version=V,code_commit=commit,relation_family_planned=len(edges),RNA_gene_family_planned=len(genes),minimum_n=8,
  input_relations_sha256=sha(edgesp),unit='unique author case within tissue; same-CAMP postselection; no genotypes',
  association_primary='Spearman rho on unchanged author processed tumor values;9999 global independent-case permutations plus1;4000 whole-case bootstrap;post-preprocessing exchangeability assumed',
  association_available='Same full relation family;require finite value in author data sheet;not certified detection mask',
  association_batch='Partial Spearman by residualizing average ranks on author batch indicators;9999 within-batch permutations;4000 batch-stratified case bootstrap',
  RNA_primary='Paired t on continuous RMA author BatchAdj matrix;43 author-case pairs;analytic t95CI and4000 bootstrap mean sensitivity;not raw counts',
  RNA_CAPT='Repeat full gene family in36 CAPT-concordant pairs;independent BH',
  missing_gene='NOT_EVALUABLE;no synonym inference or best probe selection',transform='none;reuse author processed scales',
  BH='Separate all evaluable tests per five fixed families;missing P and q remain NA',master_seed=20260922,seed='SHA256(version|20260922|family|stable_relation_or_gene_id) first8hex',
  software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__))
 js(pub/'analysis_spec.json',spec)
 rows=[]
 for i,r in edges.iterrows():
  for f in ['REL_TUMOR_PRIMARY','REL_TUMOR_AVAILABLE','REL_TUMOR_BATCH']:
   d=base(out,f,r.gene);d.update(feature_id=r.feature_id,relation_id=r.relation_id,stable_gene_id=r.stable_gene_id,metabolite_key=r.metabolite_key,metabolite_name=r.metabolite_name,n_tumor_linked=len(m))
   if r.gene not in rna.index:d.update(reason='exact_HGNC_gene_symbol_not_measured')
   else:
    x=matrix.loc[r.metabolite_name,m.MetabID].to_numpy(float);y=rna.loc[r.gene,m.RNAID].to_numpy(float);keep=np.isfinite(x)&np.isfinite(y)
    if f=='REL_TUMOR_AVAILABLE':keep &= np.isfinite(raw.loc[r.metabolite_name,m.MetabID].to_numpy(float))
    d.update(assoc(x[keep],y[keep],m.Identifier.to_numpy()[keep],f,r.relation_id,partial=f=='REL_TUMOR_BATCH'))
   rows.append(d)
  if (i+1)%20==0:print('relations_completed',i+1,flush=True)
 assocd=bh(pd.DataFrame(rows));save(assocd,pub/'tumor_association_all_families.tsv')
 rn=[]
 for _,r in genes.iterrows():
  for f in ['RNA_PAIRED_CURRENT','RNA_PAIRED_CAPT']:
   d=base(out,f,r.gene);d.update(stable_gene_id=r.stable_gene_id,expression_scale='author_RMA_BatchAdj_continuous',model_formula='RNA_T-RNA_N ~ intercept',effect_type='mean_paired_RNA_difference_author_scale',ci_method='paired_t_analytic95',ci_estimand='mean_paired_RNA_difference',p_method='two_sided_paired_t',seed=seed(f,r.stable_gene_id),bootstrap_valid=0)
   z=pairs if f=='RNA_PAIRED_CURRENT' else pairs[pairs.CAPT_concordant.eq('True')]
   if r.gene not in rna.index:d.update(reason='exact_HGNC_gene_symbol_not_measured')
   else:
    x=rna.loc[r.gene,z.RNAID_tumor].to_numpy(float);y=rna.loc[r.gene,z.RNAID_normal].to_numpy(float);keep=np.isfinite(x)&np.isfinite(y);delta=x[keep]-y[keep];n=len(delta)
    d.update(n=n,n_reference=n,n_pairs_total=len(z),n_up=int((delta>0).sum()),n_down=int((delta<0).sum()),n_equal=int((delta==0).sum()),up_fraction=float((delta>0).mean()) if n else np.nan,mean_delta=float(delta.mean()) if n else np.nan,median_delta=float(np.median(delta)) if n else np.nan)
    if n<8:d['reason']='fewer_than8_complete_pairs'
    elif np.std(delta,ddof=1)==0:d['reason']='zero_paired_difference_variance'
    else:
     eff=float(delta.mean());se=float(np.std(delta,ddof=1)/np.sqrt(n));ts=eff/se;p=float(2*stats.t.sf(abs(ts),n-1));assert abs(p-stats.ttest_rel(x[keep],y[keep]).pvalue)<1e-12
     ci=stats.t.ppf(.975,n-1)*se;rng=np.random.default_rng(d['seed']);bs=delta[rng.integers(n,size=(4000,n))].mean(1)
     d.update(effect=eff,ci_lower=eff-ci,ci_upper=eff+ci,p_value=p,t_statistic=ts,degrees_freedom=n-1,bootstrap_mean_lower=float(np.quantile(bs,.025)),bootstrap_mean_upper=float(np.quantile(bs,.975)),bootstrap_valid=4000,status='DONE',reason='NA')
   rn.append(d)
 rnad=bh(pd.DataFrame(rn));save(rnad,pub/'paired_RNA_all_families.tsv');save(rnad[rnad.test_family.eq('RNA_PAIRED_CURRENT')&rnad.p_value.lt(.05)],pub/'RNA_P005_view.tsv')
 coverage=genes.copy();coverage['RNA_exact_symbol_measured']=coverage.gene.isin(rna.index);save(coverage,pub/'RNA_identity_coverage.tsv')
 summary={f:dict(planned=len(d),evaluable=int(d.p_value.notna().sum()),p_lt005=int(d.p_value.lt(.05).sum()),q_lt005=int(d.q_value.lt(.05).sum())) for f,d in pd.concat([assocd,rnad],ignore_index=True).groupby('test_family')}
 js(pub/'validation.json',dict(status='DONE',families=summary,relation_unique=True,RNA_gene_unique=True,BH_independent_check=True,RNA_t_independent_check=True,Spearman_independent_check=True,original_values_modified=False))
 save(pd.DataFrame([dict(source_id=p.name,path_or_url=str(p),sha256=sha(p)) for p in [edgesp,mp,pp,rp,xp,Path(__file__)]]),pub/'source_manifest.tsv')
 (pub/'README_CN.md').write_text('# PRAD全关系肿瘤关联与全候选RNA背景\n\n问题：直接关系是否在PRAD肿瘤内共变，候选RNA是否在配对肿瘤中改变？\n\n输入：全部161关系、101基因、91个唯一作者肿瘤case与43对RNA；CAPT一致子集36对。\n\n实际结果：\n\n'+'\n'.join('- '+f+': '+str(v) for f,v in summary.items())+'\n\n新手解释：相关系数描述两项测量共同变化；正负相关不是合成或消耗方向。RNA变化不是酶活或通量。每个检验族单独BH，原CAMP q没有移植。\n\n限制：同队列筛选后的探索性分析，不是独立验证；作者case无基因型核实。作者批次校正后仍给批次内置换的partial Spearman敏感性，不能将批次影响隐藏。\n\n当前决定：全部基因继续单细胞来源，不按本表显著性删除；未测项保持不可评估。\n\n下一步：公开单细胞数据来源与供者标签核查；整合全部候选。\n\n复现：python prad_patient_v1.py --out AUDITED_RUN --code-commit COMMIT；读取同一服务器私有连接和固定映射输入。\n',encoding='utf-8')
 save(pd.DataFrame([dict(file=p.name,sha256=sha(p)) for p in sorted(pub.iterdir()) if p.name!='checksums.tsv']),pub/'checksums.tsv')
 print(json.dumps(summary),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--code-commit',required=True);a=ap.parse_args();main(a.out,a.code_commit)
