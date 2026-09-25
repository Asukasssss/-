"""Resolve GBM tissue conflict against pinned CPTAC loader; unpaired RNA background."""
from pathlib import Path
import argparse,requests,json,hashlib
import numpy as np,pandas as pd
from scipy import stats
import statsmodels.api as sm
from gbm_discovery_v1 import ROOT,SRC,PREFIX,sha,save,js,bh,check
V='gbm_identity_RNA_v2'
def main(out,commit):
 pub=out/'public/03_PATIENT_v2';pub.mkdir(exist_ok=False);audit=out/'public/01_CAMP_v2';audit.mkdir(exist_ok=False)
 url='https://raw.githubusercontent.com/PayneLab/cptac/v0.9.7/cptac/gbm.py';sp=out/'source/cptac_v0.9.7_gbm.py';r=requests.get(url,timeout=30);r.raise_for_status();sp.write_bytes(r.content)
 rule='sample_status_col = np.where(clinical.index.str.startswith("PT"), "Normal", "Tumor")';assert rule in sp.read_text()
 mp=out/'private/sample_identity_audit_private.tsv';m=pd.read_csv(mp,sep='\t');assert m.case_id.is_unique;m['CPTAC_loader_tissue']=np.where(m.case_id.str.startswith('PT'),'Normal','Tumor');assert m.CPTAC_loader_tissue.eq(m.TN).all();assert m.loc[m.tissue_conflict,'TN'].eq('Normal').all();save(m,out/'private/sample_identity_resolved_v2.tsv')
 summary=dict(status='DONE',mapped_author_cases=len(m),master_TN_matches_explicit_CPTAC_loader=int(m.CPTAC_loader_tissue.eq(m.TN).sum()),processed_sampleanno_conflicts=int(m.tissue_conflict.sum()),conflict_resolution='CAMP explicit master TN retained; independent CPTAC maintained package explicit PT tissue rule corroborates; derived sampleanno GROUP not used',source_loader_url=url,source_loader_sha256=sha(sp),paired_design=False,normal_reference='unmatched GTEx normal brain per Wang2021',original_files_modified=False,genotype_verification=False)
 js(audit/'validation.json',summary);js(audit/'analysis_spec.json',dict(version=V,code_commit=commit,rule=rule,decision='follow explicit CAMP tissue labels independently corroborated by maintained CPTAC loader;do not infer pairing from any ID',supersedes_interpretation='v1 normal conflict flag;v1 files preserved unchanged'))
 save(pd.DataFrame([dict(audit_item=k,value=v) for k,v in summary.items()]),audit/'sample_audit_summary.tsv')
 gp=out/'source/genes_unique.tsv';genes=pd.read_csv(gp,sep='\t');rp=SRC/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/m.RNAFile.iloc[0];rna=pd.read_csv(rp,index_col=0)
 spec=dict(version=V,code_commit=commit,genes_planned=len(genes),primary='two-sided Welch t on log2(author_TPM+1);analytic Welch mean-difference CI and4000 separate-group bootstrap',transform_reason='author TPM is nonnegative unlogged normalized abundance; log2(TPM+1) predeclared for exploratory continuous RNA background; not a raw-count differential model',effect='mean log2(TPM+1) tumor-minus-normal difference;not exact log2 fold-change',design='unpaired;74 tumor versus6 GTEx reference',minimum_group_n=5,minimum_group_reason='small normal reference;explicit exploratory unpaired adaptation;not paired n8 requirement',sensitivity='OLS log2(TPM+1) ~ tissue + age + sex;HC3 robust covariance with residual-df t;complete cases and full rank required',families=['RNA_UNPAIRED_CURRENT','RNA_UNPAIRED_AGE_SEX'],BH='separate all evaluable142 candidates in each family',normal_conflict_resolved_by=str(sp.name),master_seed=20260925,bootstrap=4000,paired_RNA='NOT_EVALUABLE_no_matched_normal_design',limitation='cohort/source and tissue confounded;not formal validation;TPM is compositional;current candidate selection from same cohort')
 js(pub/'analysis_spec.json',spec);rows=[];tum=m.TN.eq('Tumor').to_numpy();age=pd.to_numeric(m.age,errors='coerce').to_numpy();sex=m.gender.map({'Female':0,'Male':1}).to_numpy(float)
 for _,g in genes.iterrows():
  for fam in spec['families']:
   b=dict.fromkeys(PREFIX,np.nan);b.update(cancer='GBM',cohort='CAMP_GBM',stage_id='03_PATIENT',run_id=out.name,analysis_version=V,analysis_type=fam,gene=g.gene,stable_gene_id=g.stable_gene_id,unit='unpaired_unique_author_case',effect_type='mean_log2TPM1_difference',test_family=fam,status='NOT_EVALUABLE',reason='exact_symbol_missing',source_id='CAMP_GBM_TPM',same_cohort_postselection=True,normal_source_confounding=True)
   if g.gene in rna.index:
    v=rna.loc[g.gene,m.RNAID].to_numpy(float);assert np.nanmin(v)>=0;z=np.log2(v+1);keep=np.isfinite(z)
    if fam.endswith('AGE_SEX'):keep &= np.isfinite(age)&np.isfinite(sex)
    x=z[keep&tum];y=z[keep&~tum];b.update(n=len(x),n_reference=len(y))
    if min(len(x),len(y))<5:b['reason']='fewer_than5_in_group'
    elif fam=='RNA_UNPAIRED_CURRENT':
     ax=x.var(ddof=1)/len(x);ay=y.var(ddof=1)/len(y);se=np.sqrt(ax+ay)
     if se>0:
      df=(ax+ay)**2/(ax**2/(len(x)-1)+ay**2/(len(y)-1));eff=x.mean()-y.mean();p=float(2*stats.t.sf(abs(eff/se),df));assert abs(p-stats.ttest_ind(x,y,equal_var=False).pvalue)<1e-11;ci=stats.t.ppf(.975,df)*se;seed=int(hashlib.sha256((V+'|'+g.stable_gene_id).encode()).hexdigest()[:8],16);rng=np.random.default_rng(seed);boot=x[rng.integers(len(x),size=(4000,len(x)))].mean(1)-y[rng.integers(len(y),size=(4000,len(y)))].mean(1)
      b.update(effect=float(eff),ci_lower=float(eff-ci),ci_upper=float(eff+ci),p_value=p,status='DONE',reason='NA',degrees_freedom=df,p_method='Welch_t',ci_method='analytic_Welch_t95',bootstrap_mean_lower=float(np.quantile(boot,.025)),bootstrap_mean_upper=float(np.quantile(boot,.975)),bootstrap_valid=4000,seed=seed)
     else:b['reason']='zero_total_variance'
    else:
     C=np.column_stack([np.ones(keep.sum()),tum[keep].astype(int),age[keep],sex[keep]])
     if np.linalg.matrix_rank(C)==4:
      fit=sm.OLS(z[keep],C).fit(cov_type='HC3',use_t=True);ci=fit.conf_int()[1];b.update(effect=float(fit.params[1]),ci_lower=float(ci[0]),ci_upper=float(ci[1]),p_value=float(fit.pvalues[1]),status='DONE',reason='NA',effect_type='age_sex_adjusted_log2TPM1_difference',p_method='OLS_HC3_t',ci_method='HC3_t95',degrees_freedom=float(fit.df_resid))
     else:b['reason']='age_sex_design_not_full_rank'
   rows.append(b)
 d=bh(pd.DataFrame(rows));save(d,pub/'RNA_unpaired_all_families.tsv');save(d[d.test_family.eq('RNA_UNPAIRED_CURRENT')],pub/'RNA_unpaired.tsv');save(d[d.test_family.eq('RNA_UNPAIRED_CURRENT')&d.p_value.lt(.05)],pub/'RNA_P005_view.tsv')
 summaryRNA={f:dict(planned=len(z),evaluable=int(z.p_value.notna().sum()),P005=int(z.p_value.lt(.05).sum()),q005=int(z.q_value.lt(.05).sum())) for f,z in d.groupby('test_family')};js(pub/'validation.json',dict(status='DONE',families=summaryRNA,BH_independent_check=True,Welch_independent_check=True,normal_identity_source_resolved=True,paired_RNA='NOT_EVALUABLE',original_statistics_modified=False))
 man=pd.DataFrame([dict(source_id=p.name,path_or_url=url if p==sp else str(p),sha256=sha(p)) for p in [sp,mp,gp,rp,Path(__file__)]]);save(man,audit/'source_manifest.tsv');save(man,pub/'source_manifest.tsv')
 text='# GBM正常身份修订与非配对RNA背景\n\n## 本轮问题\n解决正常标签冲突并按真实非配对设计补充RNA背景。\n\n## 输入与范围\n'+json.dumps(summary,ensure_ascii=False,indent=2)+'\n\n## 实际结果\n'+json.dumps(summaryRNA,indent=2)+'\n\n## 新手解释\nCPTAC维护包对PT标本有明确正常定义；本轮沿用已有作者MasterMapping，记录processed sampleanno页错误，不更改源文件。不按编号猜配对。RNA效应是log2(TPM+1)均值差，不称酶活或精确表达倍数。\n\n## 限制/反证\n正常仅6例、GTEx与肿瘤来源混杂。Welch/年龄性别HC3分析仅为探索背景，不能消除来源混杂，不是原始counts负二项DE。无配对RNA设计。\n\n## 当前决定\n上游条件性身份标记在本修订审计中已解决；此前文件保持历史记录，整合表引用本次结论。\n\n## 下一步\n完成全基因两研究细胞来源并整合。\n\n## 复现\npython gbm_identity_RNA_v2.py --out AUDITED_RUN --code-commit COMMIT。\n'
 for p in [audit,pub]:(p/'README_CN.md').write_text(text,encoding='utf-8');check(p)
 print(json.dumps(summaryRNA),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--code-commit',required=True);a=p.parse_args();main(a.out,a.code_commit)
