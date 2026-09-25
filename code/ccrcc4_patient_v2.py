"""Patient-region means, full relationship families, explicit treatment sensitivity."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
import argparse,json,hashlib,platform
from pathlib import Path
import numpy as np,pandas as pd
from scipy import stats
from ccrcc_stats_v1 import PREFIX,save,js,sha,checksums,adjust
R=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716');S=R/'data/candidates/camp_primary_tissue_multicancer';V='ccrcc4_patient_v2'
def seed(key):return int(hashlib.sha256((V+'|20260925|'+key).encode()).hexdigest()[:8],16)
def rank_rows(x):
 order=np.argsort(x,axis=1,kind='stable');v=np.take_along_axis(x,order,axis=1);n=x.shape[1];positions=np.arange(1,n+1)[None,:]
 starts=np.concatenate([np.ones((len(x),1),bool),v[:,1:]!=v[:,:-1]],axis=1)
 ends=np.concatenate([v[:,:-1]!=v[:,1:],np.ones((len(x),1),bool)],axis=1)
 left=np.maximum.accumulate(np.where(starts,positions,0),axis=1)
 right=np.minimum.accumulate(np.where(ends,positions,n+1)[:,::-1],axis=1)[:,::-1]
 ranks=np.empty(x.shape,float);np.put_along_axis(ranks,order,(left+right)/2,axis=1)
 assert np.array_equal(ranks[:8],stats.rankdata(x[:8],axis=1))
 return ranks
def row(out,fam,gene,key='NA',name='NA'):
 d=dict.fromkeys(PREFIX,np.nan);d.update(cancer='ccRCC',cohort='ccRCC4',stage_id='03_PATIENT',run_id=out.name,analysis_version=V,analysis_type=fam,metabolite_key=key,metabolite_name=name,gene=gene,unit='author_patient',test_family=fam,status='NOT_EVALUABLE',reason='gene_not_measured',source_id='CAMP_ccRCC4_SUBJECT_ID');return d
def corr(x,y,key,strata=None):
 n=len(x);d=dict(n=n,n_reference=n,effect_type='Spearman_rho',B=9999,bootstrap_total=4000,bootstrap_valid=0,seed=seed(key),p_method='9999_patient_permutations_plus1' if strata is None else '9999_within_exact_author_treatment_permutations_plus1')
 if n<8 or len(np.unique(x))<2 or len(np.unique(y))<2:
  d.update(status='NOT_EVALUABLE',reason='fewer_than8_patients_or_constant_variable');return d
 rng=np.random.default_rng(d['seed']);rx=stats.rankdata(x);ry=stats.rankdata(y);a=rx-rx.mean();b=ry-ry.mean();den=np.sqrt(np.sum(a*a)*np.sum(b*b));rho=np.dot(a,b)/den
 assert abs(rho-stats.spearmanr(x,y).statistic)<1e-12
 perms=np.tile(np.arange(n),(9999,1))
 groups=[np.arange(n)] if strata is None else [np.flatnonzero(strata==z) for z in np.unique(strata)]
 for ix in groups:perms[:,ix]=ix[np.argsort(rng.random((9999,len(ix))),axis=1)]
 ps=(b[perms]@a)/den;p=(1+int((np.abs(ps)>=abs(rho)-1e-12).sum()))/10000
 boot=np.empty((4000,n),dtype=int)
 for ix in groups:boot[:,ix]=rng.choice(ix,size=(4000,len(ix)),replace=True)
 xx=rank_rows(x[boot]);yy=rank_rows(y[boot]);xx-=xx.mean(1,keepdims=True);yy-=yy.mean(1,keepdims=True)
 bd=np.sqrt((xx*xx).sum(1)*(yy*yy).sum(1));valid=bd>0;bs=(xx[valid]*yy[valid]).sum(1)/bd[valid]
 d.update(effect=float(rho),rho=float(rho),p_value=p,status='DONE',reason='NA',bootstrap_valid=int(valid.sum()),ci_method='pointwise_whole_patient_percentile_bootstrap' if strata is None else 'pointwise_treatment_stratified_whole_patient_percentile_bootstrap')
 if len(bs)>=3600:d.update(ci_lower=float(np.quantile(bs,.025)),ci_upper=float(np.quantile(bs,.975)))
 else:d['reason']='bootstrap_valid_below90pct_CI_not_evaluable'
 return d
def main(out,commit):
 with (out/'.patient_running').open('x') as f:f.write(V)
 pub=out/'public/03_PATIENT';pub.mkdir(exist_ok=False)
 mp=out/'private/sample_identity_audit_private.tsv';m=pd.read_csv(mp,sep='\t',dtype=str)
 ep=out/'source/direct_relations.tsv';gp=out/'source/gene_membership.tsv';edges=pd.read_csv(ep,sep='\t');genes=pd.read_csv(gp,sep='\t')
 rp=S/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/m.RNAFile.iloc[0];rna=pd.read_csv(rp,index_col=0)
 assert rna.index.is_unique and rna.columns.is_unique and m.RNAID.isin(rna.columns).all()
 assert np.isfinite(rna.loc[:,m.RNAID].to_numpy()).all() and rna.loc[:,m.RNAID].ge(0).all().all()
 metp=S/'processed_metabolomics/PreprocessedData_ccRCC4.xlsx';mat=pd.read_excel(metp,sheet_name='data_imputed',index_col=0);raw=pd.read_excel(metp,sheet_name='data',index_col=0)
 for z in [mat,raw]:z.index=z.index.astype(str).str.strip()
 spec=dict(version=V,code_commit=commit,unit='explicit author SUBJECT_ID, regions averaged separately within tissue',
  RNA_input='author TPM gene-symbol file;nonnegative continuous abundance,not raw counts',RNA_transform='log2(1+TPM) per specimen before patient-region mean;new explicit analysis transform',
  metabolite_transform='none;author processed values',families=['REL_TUMOR_PRIMARY','REL_TUMOR_AVAILABLE','REL_TUMOR_TREATMENT_STRATIFIED','RNA_PAIRED_CURRENT','RNA_PAIRED_HISTORY_ONLY'],
  association='Spearman of patient tumor mean processed metabolite versus patient mean log2(1+TPM)',
  p_method='9999 two-sided patient permutations plus1;prespecified sensitivity permutes within exact author TREATMENT',
  CI='4000 whole-patient bootstrap, percentile;stratified sensitivity resamples within exact treatment',
  RNA_test='paired t test on patient region means; t CI plus4000 patient bootstrap;zero variance or n<8 NOT_EVALUABLE',
  BH='all evaluable tests independently per family, including all current relations/current genes;history-only gene family separate',
  availability='all regions included in patient tumor mean must have author data values;same region mixture retained',
  selection='same_cohort_postselection,not independent validation',master_seed=20260925,software=dict(python=platform.python_version(),pandas=pd.__version__,numpy=np.__version__))
 js(pub/'analysis_spec.json',spec)
 tumor=m[m.TN.eq('Tumor')];groups=tumor.groupby('SUBJECT_ID');pats=sorted(groups.groups)
 treatment=groups.TREATMENT.first().reindex(pats).to_numpy();metids={k:z.MetabID.tolist() for k,z in groups};rnaids={k:z.RNAID.tolist() for k,z in groups}
 pairs=sorted(set(pats)&set(m[m.TN.eq('Normal')].SUBJECT_ID));ng=m[m.TN.eq('Normal')].groupby('SUBJECT_ID').RNAID.apply(list).to_dict()
 rows=[]
 for _,e in edges.iterrows():
  for fam in ['REL_TUMOR_PRIMARY','REL_TUMOR_AVAILABLE','REL_TUMOR_TREATMENT_STRATIFIED']:
   d=row(out,fam,e.gene,e.metabolite_key,e.metabolite_name);d.update(feature_id=e.feature_id,relation_id=e.relation_id,stable_gene_id=e.stable_gene_id,n_tumor_linked=len(pats),same_cohort_postselection=True)
   if e.gene in rna.index:
    x=np.array([mat.loc[e.metabolite_name,metids[p]].mean() for p in pats]);y=np.array([np.log2(1+rna.loc[e.gene,rnaids[p]].to_numpy(float)).mean() for p in pats])
    mask=np.isfinite(x)&np.isfinite(y)
    if fam=='REL_TUMOR_AVAILABLE':mask&=np.array([np.isfinite(raw.loc[e.metabolite_name,metids[p]]).all() for p in pats])
    d.update(corr(x[mask],y[mask],fam+'|'+e.relation_id,treatment[mask] if fam.endswith('STRATIFIED') else None))
   rows.append(d)
  if len(rows)%75==0:print('RELATIONS',len(rows)//3,flush=True)
 associations=adjust(pd.DataFrame(rows));save(associations,pub/'association_results.tsv')
 rows=[]
 for _,g in genes.iterrows():
  fam='RNA_PAIRED_CURRENT' if g.current_pool else 'RNA_PAIRED_HISTORY_ONLY';d=row(out,fam,g.gene);d.update(stable_gene_id=g.stable_gene_id,current_pool=bool(g.current_pool),history_only=bool(g.history_only),n_pairs_total=len(pairs),expression_scale='patient_mean_log2_1plus_author_TPM',effect_type='mean_paired_difference_log2_1plus_TPM')
  if g.gene in rna.index:
   delta=np.array([np.log2(1+rna.loc[g.gene,rnaids[p]].to_numpy(float)).mean()-np.log2(1+rna.loc[g.gene,ng[p]].to_numpy(float)).mean() for p in pairs]);delta=delta[np.isfinite(delta)];n=len(delta)
   d.update(n=n,n_reference=n,n_up=int((delta>0).sum()),n_down=int((delta<0).sum()),n_equal=int((delta==0).sum()),effect=float(delta.mean()),mean_delta=float(delta.mean()),median_delta=float(np.median(delta)),up_fraction=float((delta>0).mean()),seed=seed(fam+'|'+g.gene))
   if n>=8 and delta.std(ddof=1)>0:
    se=stats.sem(delta);tval=delta.mean()/se;p=2*stats.t.sf(abs(tval),n-1);assert abs(p-stats.ttest_1samp(delta,0).pvalue)<1e-12
    ci=stats.t.interval(.95,n-1,loc=delta.mean(),scale=se);rng=np.random.default_rng(d['seed']);bs=delta[rng.integers(n,size=(4000,n))].mean(1)
    d.update(p_value=float(p),ci_lower=ci[0],ci_upper=ci[1],p_method='paired_t_on_patient_means',ci_method='t_df_nminus1_pointwise',bootstrap_mean_lower=np.quantile(bs,.025),bootstrap_mean_upper=np.quantile(bs,.975),status='DONE',reason='NA')
   else:d.update(reason='fewer_than8_pairs_or_zero_delta_variance')
  rows.append(d)
 rr=adjust(pd.DataFrame(rows));save(rr,pub/'rna_results.tsv')
 assert len(associations)==3*len(edges) and len(rr)==len(genes) and not associations.duplicated(['relation_id','test_family']).any()
 summary={}
 for name,d in [('association',associations),('RNA',rr)]:summary[name]=d.groupby('test_family').apply(lambda z:dict(planned=len(z),evaluable=int(z.p_value.notna().sum()),p005=int(z.p_value.lt(.05).sum()),q005=int(z.q_value.lt(.05).sum()))).to_dict()
 val=dict(status='DONE',tumor_patients=len(pats),paired_patients=len(pairs),families=summary,independent_rho_and_t_checks=True,BH_formula_independently_checked=True,all_candidates_retained=True)
 js(pub/'validation.json',val);save(pd.DataFrame([dict(source_id=p.name,path_or_url=str(p),sha256=sha(p),access_scope='PRIVATE' if p in [mp,rp,metp] else 'CODE_OR_SUMMARY') for p in [mp,ep,gp,rp,metp,Path(__file__)] ]),pub/'source_manifest.tsv')
 (pub/'README_CN.md').write_text(f'''# ccRCC4全候选患者关联与RNA背景\n\n问题：肿瘤患者内代谢物和对应基因是否协变？配对RNA背景如何？\n\n输入范围：{len(pats)}位肿瘤患者、{len(pairs)}对；{len(edges)}条当前关系、{len(genes)}个当前/历史并集基因。多区域先取患者内均值；RNA明确采用log2(1+作者TPM)，不是原始计数差异分析。\n\n实际结果：见validation.json各检验族计数，完整结果保留不可评估项。\n\n新手解释：rho是秩相关，不是倍数；RNA差异和相关检验独立，q分别校正。\n\n限制：12位患者的小样本和治疗混合使结果不精确；主分析未调整治疗，另列在作者精确治疗类别内置换的敏感性，同患者均值无法刻画区域内机制。候选来自同队列筛选，不是独立验证。\n\n当前决定：全候选进入单细胞来源，不按显著性删基因。\n\n下一步：表达来源与整合。\n\n复现：python ccrcc4_patient_v2.py --out SERVER_RUN --code-commit COMMIT。\n''',encoding='utf-8')
 checksums(pub);(out/'.patient_running').rename(out/'.patient_done');print(json.dumps(val),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--code-commit',required=True);a=p.parse_args();main(a.out,a.code_commit)


