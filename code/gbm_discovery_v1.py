"""GBM server-only discovery: exact author mapping, frozen effects, conditional normal labels."""
from pathlib import Path
import argparse,hashlib,json,platform
import numpy as np,pandas as pd,scipy
from scipy import stats
from statsmodels.stats.multitest import multipletests
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
SRC=ROOT/'data/candidates/camp_primary_tissue_multicancer'
V='gbm_discovery_v1'
PREFIX='cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA')
def js(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
def seed(f,k):return int(hashlib.sha256((V+'|20260925|'+f+'|'+k).encode()).hexdigest()[:8],16)
def check(out):save(pd.DataFrame([dict(file=p.name,sha256=sha(p)) for p in sorted(out.iterdir()) if p.is_file() and p.name!='checksums.tsv']),out/'checksums.tsv')
def bh(d):
 for f,ix in d.groupby('test_family').groups.items():
  ok=d.index.isin(ix)&d.p_value.notna();p=d.loc[ok,'p_value'].to_numpy(float);d.loc[ix,'family_n_planned']=len(ix);d.loc[ix,'family_n_evaluable']=len(p)
  if len(p):
   q=multipletests(p,method='fdr_bh')[1];order=np.argsort(p);ref=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);assert np.allclose(q[order],ref,atol=1e-13);d.loc[ok,'q_value']=q
 return d

def main(out,commit):
 out.resolve().relative_to(ROOT/'results/collaborative/GBM/A')
 with (out/'.running').open('x') as f:f.write(V)
 for x in ['private','public/01_CAMP','public/04_ROBUSTNESS']:(out/x).mkdir(parents=True,exist_ok=True)
 pub=out/'public/01_CAMP';rob=out/'public/04_ROBUSTNESS'
 mp=SRC/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv';xp=SRC/'processed_metabolomics/PreprocessedData_GBM.xlsx';ap=out/'source/GBM_original.xlsx'
 m=pd.read_csv(mp,dtype=str);m=m[m.Dataset.eq('GBM')].copy();a=pd.read_excel(ap,sheet_name='sampleanno');assert a.case_id.is_unique
 sa=pd.read_excel(xp,sheet_name='sampleanno',index_col=0);m['processed_annotation']=m.MetabID.map(sa.GROUP);m['tissue_conflict']=m.TN.str.upper().ne(m.processed_annotation)
 assert all(m[c].is_unique and m[c].notna().all() for c in ['RNAID','MetabID','CommonID']);assert m.RNAID.eq(m.MetabID).all();assert m.MetabID.isin(a.case_id).all()
 m=m.merge(a,left_on='MetabID',right_on='case_id',validate='one_to_one');m['unit_verification']='exact_author_case_id;no_genotype_check';m['link_level']='author_case_matched_omics;aliquot_not_separately_verified'
 save(m,out/'private/sample_identity_audit_private.tsv');t=m[m.TN.eq('Tumor')&~m.tissue_conflict].copy();save(t,out/'private/tumor_multiomics_map_private.tsv')
 rp=SRC/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/m.RNAFile.iloc[0];rna=pd.read_csv(rp,index_col=0);assert m.RNAID.isin(rna.columns).all() and rna.index.is_unique
 raw=pd.read_excel(xp,sheet_name='data',index_col=0);mat=pd.read_excel(xp,sheet_name='data_imputed',index_col=0);raw.index=raw.index.str.strip();mat.index=mat.index.str.strip();ann=pd.read_excel(xp,sheet_name='metanno');assert mat.index.is_unique and raw.index.equals(mat.index)
 oldp=ROOT/'results/tables/cohort_effects.tsv';old=pd.read_csv(oldp,sep='\t');old=old[old.dataset.eq('GBM')].copy();assert old.feature_name.is_unique
 cp=ROOT/'results/tables/cancer_effects.tsv';frozen=pd.read_csv(cp,sep='\t');frozen=frozen[frozen.cancer.eq('GBM')].copy()
 save(old,pub/'historical_cohort_effects_preserved.tsv');save(frozen,pub/'historical_cancer_effects_preserved.tsv')
 retained=[]
 for typ in ['Tumor','Normal']:
  d=pd.read_excel(xp,sheet_name='metabo_imputed_filtered_'+typ,index_col=0);d.index=d.index.str.strip();assert set(d.columns)==set(m.loc[m.TN.eq(typ),'MetabID']);assert np.allclose(d.values,mat.loc[d.index,d.columns].values);retained.append(set(d.index))
 assert set(old.feature_name)==retained[0]&retained[1]
 inv=pd.DataFrame({'metabolite_name':mat.index});inv['feature_id']=inv.metabolite_name.map(lambda z:'GBM_FEATURE:'+hashlib.sha256(z.encode()).hexdigest()[:16]);inv['retained_both_tissues']=inv.metabolite_name.isin(old.feature_name)
 inv=inv.merge(old[['feature_name','metabolite_key','identity_confidence']],left_on='metabolite_name',right_on='feature_name',how='left',validate='one_to_one');inv['status']=np.where(inv.retained_both_tissues,'NEEDS_REVIEW','NOT_EVALUABLE');inv['reason']=np.where(inv.retained_both_tissues,'six_normal_tissue_labels_conflict_between_author_files','not_retained_in_both_tissue_filters');save(inv,pub/'feature_inventory.tsv')
 spec=dict(analysis_version=V,code_commit=commit,primary='reuse exact historical Mann-Whitney P/effect on same inputs; no frozen statistic overwritten; conditional on master TN',workpool='historical full-cohort P<0.05; exploratory and conditional normal-label conflict',normal_label_conflicts=int(m.tissue_conflict.sum()),normal_discovery_status='NEEDS_REVIEW',patient_unit='unique explicit author case_id',paired_analysis='NOT_EVALUABLE;unmatched normal brain reference',metabolite_scale='original README log2 transformed globally median normalized; no new transform',minimum_sensitivity_group_n=5,minimum_n_reason='unpaired adaptation fixed before new sensitivity statistics; only six reference normals; paired/association n8 rule not reused for group size',sensitivity='Mann-Whitney asymptotic,two-sided,tie and continuity correction;author data finite-value mask;separate full family BH;4000 within-group mean bootstrap',bootstrap=4000,master_seed=20260925,mask_semantics='author_available_value_mask;not verified original detection mask',clinical_sensitivity='exclude author pathology FAIL from tumor associations;age and sex adjustment when estimable',software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__),original_statistics_modified=False)
 js(pub/'analysis_spec.json',spec);js(rob/'analysis_spec.json',spec)
 rows=[];sens=[]
 for _,r in inv.iterrows():
  b=dict.fromkeys(PREFIX,np.nan);b.update(cancer='GBM',cohort='CAMP_GBM',stage_id='01_CAMP',run_id=out.name,analysis_version=V,analysis_type='historical_unpaired_reuse_conditional',metabolite_key=r.metabolite_key,metabolite_name=r.metabolite_name,feature_id=r.feature_id,gene='NA',unit='author_case',test_family='METAB_UNPAIRED_HISTORICAL',source_id='CAMP_v0.3.4_GBM',status=r.status,reason=r.reason,normal_label_conflict=True)
  if r.retained_both_tissues:
   z=old.set_index('feature_name').loc[r.metabolite_name];b.update(n=z.n_tumor,n_reference=z.n_normal,effect_type='Hedges_g_preserved',effect=z.hedges_g,ci_lower=z.hedges_g_ci_lower,ci_upper=z.hedges_g_ci_upper,p_value=z.wilcoxon_p,original_cohort_q=z.wilcoxon_fdr,statistics_reused=True)
   x=mat.loc[r.metabolite_name,m.loc[m.TN.eq('Tumor'),'MetabID']].to_numpy(float);y=mat.loc[r.metabolite_name,m.loc[m.TN.eq('Normal'),'MetabID']].to_numpy(float)
   p=stats.mannwhitneyu(x,y,alternative='two-sided',method='asymptotic',use_continuity=True).pvalue;assert abs(p-z.wilcoxon_p)<1e-10,(r.metabolite_name,p,z.wilcoxon_p)
   q=b.copy();q.update(stage_id='04_ROBUSTNESS',analysis_type='available_unpaired_conditional',test_family='METAB_UNPAIRED_AVAILABLE',effect_type='mean_difference_author_log2_scale',effect=np.nan,ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,statistics_reused=False)
   xt=raw.loc[r.metabolite_name,m.loc[m.TN.eq('Tumor'),'MetabID']].notna().to_numpy();yn=raw.loc[r.metabolite_name,m.loc[m.TN.eq('Normal'),'MetabID']].notna().to_numpy();x=x[xt];y=y[yn];q.update(n=len(x),n_reference=len(y),mask_semantics=spec['mask_semantics'])
   if min(len(x),len(y))>=5:
    rng=np.random.default_rng(seed(q['test_family'],r.feature_id));bs=x[rng.integers(len(x),size=(4000,len(x)))].mean(1)-y[rng.integers(len(y),size=(4000,len(y)))].mean(1);q.update(effect=float(x.mean()-y.mean()),ci_lower=float(np.quantile(bs,.025)),ci_upper=float(np.quantile(bs,.975)),p_value=float(stats.mannwhitneyu(x,y,alternative='two-sided',method='asymptotic',use_continuity=True).pvalue),bootstrap_valid=4000,seed=seed(q['test_family'],r.feature_id))
   else:q.update(status='NOT_EVALUABLE',reason='fewer_than5_available_in_a_group;normal_labels_also_conflicted')
   sens.append(q)
  rows.append(b)
 d=bh(pd.DataFrame(rows));ss=bh(pd.DataFrame(sens));save(d,pub/'metabolite_unpaired_all.tsv');save(d[d.p_value.lt(.05)],pub/'metabolite_workpool.tsv');save(ss,rob/'metabolite_available_sensitivity.tsv')
 summary=dict(status='PARTIAL',total_input_features=len(inv),retained_both_tissues=len(old),historical_cancer_features=len(frozen),workpool_p005=int(d.p_value.lt(.05).sum()),current_family_q005=int(d.q_value.lt(.05).sum()),tumor_cases=len(t),normal_cases=int(m.TN.eq('Normal').sum()),normal_label_conflicts=int(m.tissue_conflict.sum()),tumor_pathology_PASS=int(t.expert_path_review.eq('PASS').sum()),tumor_pathology_FAIL=int(t.expert_path_review.eq('FAIL').sum()),RNA_exact_ID_coverage=len(m),RNA_gene_symbols=len(rna),historical_P_independently_reproduced=True,BH_independently_checked=True,all_source_feature_rows_preserved=True,original_statistics_modified=False,genotype_check='NOT_EVALUABLE_no_genotype_data',sex_age='available_from_original_clinical;used_in_next_stage')
 js(pub/'validation.json',summary);js(rob/'validation.json',dict(status='PARTIAL',tests_evaluable=int(ss.p_value.notna().sum()),normal_identity_conflict_unresolved=True,BH_independently_checked=True))
 save(pd.DataFrame([dict(audit_item=k,value=v) for k,v in summary.items()]),pub/'sample_audit_summary.tsv')
 sources=[mp,xp,ap,rp,oldp,cp,Path(__file__)];manifest=pd.DataFrame([dict(source_id=p.name,path_or_url=str(p),sha256=sha(p),access_scope='SERVER_PRIVATE') for p in sources]);save(manifest,pub/'source_manifest.tsv');save(manifest,rob/'source_manifest.tsv')
 text='# GBM发现与身份核查\n\n## 本轮问题\n按BRCA主线建立GBM全量发现与候选工作池。\n\n## 输入与范围\nCAMP作者处理矩阵、明确MasterMapping、原始GBM.xlsx病例资料；仅服务器读取。\n\n## 实际结果\n'+json.dumps(summary,ensure_ascii=False,indent=2)+'\n\n## 新手解释\n原始P/效应精确复用，当前完整特征族的BH另列。不同脂质峰即使同KEGG仍保留，不合并。P<0.05工作池是探索入口。\n\n## 限制/反证\n74肿瘤与6正常来自非配对设计。6正常的MasterMapping和processed sampleanno标签冲突，尚未以独立明确组织标签解决。原临床的病例ID均能精确连接，但缺少肿瘤字段不能单独证明正常身份。因此所有依赖正常组的结果为NEEDS_REVIEW，不升级确认性证据。原研究正常来自GTEx，与肿瘤来源有混杂；缺失掩码仅为作者可用值。\n\n## 当前决定\n保留冻结结果和条件性全量工作池；不伪造配对、不自动更改组织标签。74个组织一致的肿瘤可做内部关系，另给68个病理PASS子集敏感性。\n\n## 下一步\n直接映射、肿瘤内关联、全候选细胞来源；正常依赖分析待身份解决。\n\n## 复现\npython gbm_discovery_v1.py --out NEW_SERVER_RUN --code-commit COMMIT。source放GBM_original.xlsx。\n'
 (pub/'README_CN.md').write_text(text,encoding='utf-8');(rob/'README_CN.md').write_text(text+'\n可用值敏感性按每组至少5个预定规则执行；不代表解决组织冲突。\n',encoding='utf-8');check(pub);check(rob);print(json.dumps(summary),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--code-commit',required=True);a=p.parse_args();main(a.out,a.code_commit)

