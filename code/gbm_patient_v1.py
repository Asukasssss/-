"""Server-only GBM tumor associations with explicit clinical QC and covariate sensitivity."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
from pathlib import Path
import argparse,json,hashlib,platform
import numpy as np,pandas as pd,scipy
from scipy import stats
from gbm_discovery_v1 import ROOT,SRC,PREFIX,sha,save,js,bh,check
V='gbm_patient_v1'
def seed(f,k):return int(hashlib.sha256((V+'|20260925|'+f+'|'+k).encode()).hexdigest()[:8],16)
def assoc(x,y,f,key,cov=None):
 n=len(x);s=seed(f,key);d=dict(n=n,n_reference=np.nan,seed=s,bootstrap_valid=0,ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,effect=np.nan,status='NOT_EVALUABLE',reason='fewer_than8_or_constant',effect_type='Spearman_rho',p_method='9999_permutations_plus1',ci_method='4000_whole_case_percentile_bootstrap_pointwise')
 if n<8 or np.ptp(x)==0 or np.ptp(y)==0:return d
 rx=stats.rankdata(x);ry=stats.rankdata(y);rng=np.random.default_rng(s)
 if cov is not None:
  C=np.column_stack([np.ones(n),stats.rankdata(cov[:,0]),cov[:,1]]);rank=np.linalg.matrix_rank(C)
  if rank!=3 or n<=rank+2:d['reason']='age_sex_design_not_estimable';return d
  rx=rx-C@np.linalg.lstsq(C,rx,rcond=None)[0];ry=ry-C@np.linalg.lstsq(C,ry,rcond=None)[0]
 else:rx-=rx.mean();ry-=ry.mean()
 norm=np.linalg.norm(rx)*np.linalg.norm(ry)
 if norm==0:return d
 rho=float(rx@ry/norm)
 if cov is None:
  assert abs(rho-stats.spearmanr(x,y).statistic)<1e-12
  extreme=0
  for k in range(0,9999,1000):
   ix=np.array([rng.permutation(n) for _ in range(min(1000,9999-k))]);v=ry[ix]@rx/norm;extreme+=int((abs(v)>=abs(rho)-1e-12).sum())
  p=(extreme+1)/10000
 else:
  df=n-4;p=float(2*stats.t.sf(abs(rho)*np.sqrt(df/max(1e-15,1-rho*rho)),df));d.update(effect_type='partial_rank_rho_age_sex',p_method='approximate_partial_rank_t_df_n_minus4',ci_method='NOT_COMPUTED_covariate_sensitivity',model_formula='rank(metabolite) ~ rank(age)+sex;rank(TPM) ~ rank(age)+sex',degrees_freedom=df)
 if cov is None:
  ix=rng.integers(n,size=(4000,n));a=stats.rankdata(x[ix],axis=1);b=stats.rankdata(y[ix],axis=1);a-=a.mean(1,keepdims=True);b-=b.mean(1,keepdims=True);den=np.sqrt((a*a).sum(1)*(b*b).sum(1));valid=den>0;boot=(a[valid]*b[valid]).sum(1)/den[valid];d['bootstrap_valid']=int(valid.sum())
  if valid.sum()>=3600:d.update(ci_lower=float(np.quantile(boot,.025)),ci_upper=float(np.quantile(boot,.975)))
 d.update(effect=rho,rho=rho,p_value=p,status='DONE',reason='NA');return d

def main(out,commit):
 assert (out/'.running').read_text()=='gbm_discovery_v1';pub=out/'public/03_PATIENT';pub.mkdir(exist_ok=False)
 ep=out/'source/direct_relations.tsv';edges=pd.read_csv(ep,sep='\t');genes=edges[['gene','stable_gene_id']].drop_duplicates();assert genes.gene.is_unique and not edges.duplicated(['feature_id','stable_gene_id']).any()
 mp=out/'private/tumor_multiomics_map_private.tsv';m=pd.read_csv(mp,sep='\t');assert m.case_id.is_unique and not m.tissue_conflict.any()
 rp=SRC/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/m.RNAFile.iloc[0];xp=SRC/'processed_metabolomics/PreprocessedData_GBM.xlsx';rna=pd.read_csv(rp,index_col=0);mat=pd.read_excel(xp,sheet_name='data_imputed',index_col=0);raw=pd.read_excel(xp,sheet_name='data',index_col=0);mat.index=mat.index.str.strip();raw.index=raw.index.str.strip();assert mat.index.is_unique and rna.index.is_unique
 spec=dict(version=V,code_commit=commit,planned_relations=len(edges),planned_genes=len(genes),tumor_cases=len(m),normal_specimens_used=0,families=['REL_TUMOR_PRIMARY','REL_TUMOR_AVAILABLE','REL_TUMOR_PATH_PASS','REL_TUMOR_AGE_SEX'],primary='unaltered author metabolite scale versus unaltered author TPM;Spearman;9999 permutation plus1;4000 case bootstrap',pathology_sensitivity='explicit original expert_path_review PASS only',age_sex='complete age,sex;partial rank residual Pearson;approx t n-4,not exact permutation;separate BH',minimum_n=8,master_seed=20260925,seed_derivation='SHA256(version|20260925|family|relation_id) first8hex',RNA_background_status='NOT_EVALUABLE:no paired normals and unresolved normal tissue labels;TPM not tested as raw counts',same_cohort_postselection=True,normal_discovery_gate='conditional;normal labels conflict',mask_semantics='author_available_value_mask',software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__))
 js(pub/'analysis_spec.json',spec);rows=[];age=pd.to_numeric(m.age,errors='coerce').to_numpy();sex=m.gender.map({'Female':0,'Male':1}).to_numpy(float)
 for i,r in edges.iterrows():
  for fam in spec['families']:
   d=dict.fromkeys(PREFIX,np.nan);d.update(cancer='GBM',cohort='CAMP_GBM',stage_id='03_PATIENT',run_id=out.name,analysis_version=V,analysis_type=fam,metabolite_key=r.metabolite_key,metabolite_name=r.metabolite_name,gene=r.gene,stable_gene_id=r.stable_gene_id,feature_id=r.feature_id,relation_id=r.relation_id,unit='unique_explicit_author_case',test_family=fam,status='NOT_EVALUABLE',reason='exact_gene_symbol_missing',source_id='CAMP_v0.3.4_GBM',same_cohort_postselection=True,discovery_gate='CONDITIONAL_NORMAL_LABEL_CONFLICT')
   if r.gene in rna.index:
    x=mat.loc[r.metabolite_name,m.MetabID].to_numpy(float);y=rna.loc[r.gene,m.RNAID].to_numpy(float);keep=np.isfinite(x)&np.isfinite(y)
    if fam=='REL_TUMOR_AVAILABLE':keep &= raw.loc[r.metabolite_name,m.MetabID].notna().to_numpy()
    if fam=='REL_TUMOR_PATH_PASS':keep &= m.expert_path_review.eq('PASS').to_numpy()
    if fam=='REL_TUMOR_AGE_SEX':keep &= np.isfinite(age)&np.isfinite(sex)
    cov=np.column_stack([age[keep],sex[keep]]) if fam=='REL_TUMOR_AGE_SEX' else None;d.update(assoc(x[keep],y[keep],fam,r.relation_id,cov))
   rows.append(d)
  if (i+1)%20==0:print('relations',i+1,flush=True)
 d=bh(pd.DataFrame(rows));d=d[PREFIX+[c for c in d if c not in PREFIX]];save(d,pub/'tumor_association_all_families.tsv');save(d[d.test_family.eq('REL_TUMOR_PRIMARY')],pub/'tumor_association.tsv')
 rr=[]
 for _,r in genes.iterrows():
  q=dict.fromkeys(PREFIX,np.nan);q.update(cancer='GBM',cohort='CAMP_GBM',stage_id='03_PATIENT',run_id=out.name,analysis_version=V,analysis_type='RNA_BACKGROUND',gene=r.gene,stable_gene_id=r.stable_gene_id,test_family='RNA_BACKGROUND_CURRENT',status='NOT_EVALUABLE',reason='no_paired_normal_design;normal_tissue_label_conflict_unresolved;no_RNA_differential_test_performed',source_id='CAMP_v0.3.4_GBM',RNA_exact_symbol_measured=r.gene in rna.index);rr.append(q)
 save(pd.DataFrame(rr),pub/'RNA_background_status.tsv')
 summary={f:dict(planned=len(z),evaluable=int(z.p_value.notna().sum()),P005=int(z.p_value.lt(.05).sum()),q005=int(z.q_value.lt(.05).sum())) for f,z in d.groupby('test_family')};js(pub/'validation.json',dict(status='PARTIAL',families=summary,Spearman_independent_check=True,BH_independent_check=True,no_normal_values_used=True,RNA_test_not_run=True,original_statistics_modified=False))
 save(pd.DataFrame([dict(source_id=p.name,path_or_url=str(p),sha256=sha(p)) for p in [ep,mp,rp,xp,Path(__file__),Path(__file__).with_name('gbm_discovery_v1.py')]]),pub/'source_manifest.tsv')
 (pub/'README_CN.md').write_text('# GBM全直接关系肿瘤内关联\n\n## 本轮问题\n直接代谢物—基因关系在肿瘤病例中是否共变？\n\n## 输入与范围\n'+str(len(edges))+'条关系、'+str(len(genes))+'基因、74个组织标签一致的独立作者case；原处理代谢值与TPM不再变换。\n\n## 实际结果\n'+json.dumps(summary,indent=2)+'\n\n## 新手解释\n效应为秩相关，不是倍数、酶活或合成方向。全部关系各家族独立BH。\n\n## 限制/反证\n同队列选择后分析；上游候选入口仍依赖身份待核实的正常组。关联本身不使用正常标本。病理PASS敏感性排除6个低质量病例；年龄性别partial rank为近似t检验。未控制IDH/纯度等未冻结协变量，不称独立验证。RNA配对和正常差异暂不可评估，不填假P。\n\n## 当前决定\n所有142基因保留进入单细胞，不按P筛减。\n\n## 下一步\n全候选细胞来源描述和整合。\n\n## 复现\npython gbm_patient_v1.py --out AUDITED_RUN --code-commit COMMIT；同目录需要gbm_discovery_v1.py。\n',encoding='utf-8');check(pub);print(json.dumps(summary),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--code-commit',required=True);a=p.parse_args();main(a.out,a.code_commit)
