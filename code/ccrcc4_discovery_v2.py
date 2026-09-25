"""ccRCC4 patient-level discovery; ccRCC3 remains gated by author patient mapping."""
import argparse, json, platform, sys
from pathlib import Path
import pandas as pd
import numpy as np
from ccrcc_stats_v1 import sha, save, js, checksums, paired_stats, adjust, PREFIX
R=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
S=R/'data/candidates/camp_primary_tissue_multicancer'
V='ccrcc4_discovery_v2'

def main(out,commit):
 out.resolve().relative_to(R/'results/collaborative/ccRCC/B')
 with (out/'.running').open('x') as f:f.write(V)
 private=out/'private';private.mkdir(exist_ok=False)
 pub=out/'public/01_CAMP';pub.mkdir(parents=True,exist_ok=False)
 source=out/'source/ccRCC4_original.xlsx'
 mp=S/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv'
 original=pd.read_excel(source,sheet_name='OrigScale',header=None,nrows=21)
 a=original.iloc[:,13:].T;a.columns=original.iloc[:,12].tolist()
 assert a.SAMPLE_NAME.is_unique and a.SUBJECT_ID.notna().all()
 mapping=pd.read_csv(mp,dtype=str)
 m=mapping[mapping.Dataset.eq('ccRCC4')].merge(a,left_on='MetabID',right_on='SAMPLE_NAME',validate='one_to_one')
 assert len(m)==76 and m.MetabID.is_unique and m.RNAID.is_unique
 assert m.TN.map({'Tumor':'TISSUE_TUMOR','Normal':'TISSUE_NORMAL'}).eq(m.SAMPLE_DESCRIPTION).all()
 assert m.groupby('SUBJECT_ID').TREATMENT.nunique().max()==1
 meta5path=out/'source/MetaData_M5.tsv'
 meta5=pd.read_csv(meta5path,sep='\t',dtype=str)
 meta5['MetabID']=meta5.SAMPLE_NAME.str.replace('.','-',regex=False)
 m=m.merge(meta5[['MetabID','SUBJECT_ID','HISTOLOGY']],on='MetabID',suffixes=('','_study'),validate='one_to_one')
 assert len(m)==76 and m.SUBJECT_ID.eq(m.SUBJECT_ID_study).all()
 bad=set(m.loc[~m.Histology.eq('Clear Cell')|~m.HISTOLOGY.eq('Clear Cell'),'SUBJECT_ID'])
 excluded=m[m.SUBJECT_ID.isin(bad)].copy();m=m[~m.SUBJECT_ID.isin(bad)].copy()
 assert m.groupby('SUBJECT_ID').Histology.nunique().max()==1
 save(m,private/'sample_identity_audit_private.tsv');save(excluded,private/'histology_excluded_private.tsv')
 paired=sorted(set(m[m.TN.eq('Tumor')].SUBJECT_ID)&set(m[m.TN.eq('Normal')].SUBJECT_ID))
 save(pd.DataFrame({'SUBJECT_ID':paired}),private/'paired_patients_private.tsv')
 path=S/'processed_metabolomics/PreprocessedData_ccRCC4.xlsx';x=pd.ExcelFile(path)
 mat=pd.read_excel(x,sheet_name='data_imputed',index_col=0);raw=pd.read_excel(x,sheet_name='data',index_col=0)
 for z in [mat,raw]:z.index=z.index.astype(str).str.strip()
 assert mat.index.is_unique and mat.columns.is_unique and mat.index.equals(raw.index) and mat.columns.equals(raw.columns)
 assert m.MetabID.isin(mat.columns).all()
 assert np.allclose(mat.values[np.isfinite(raw.values)],raw.values[np.isfinite(raw.values)])
 ft=pd.read_excel(x,sheet_name='metabo_imputed_filtered_Tumor',index_col=0)
 fn=pd.read_excel(x,sheet_name='metabo_imputed_filtered_Normal',index_col=0)
 retained=set(ft.index.astype(str).str.strip())&set(fn.index.astype(str).str.strip())
 assert set(m.MetabID)<=set(ft.columns)|set(fn.columns)
 for tissue,d in [('Tumor',ft),('Normal',fn)]:
  d.index=d.index.astype(str).str.strip()
  assert np.allclose(d.values,mat.loc[d.index,d.columns].values,equal_nan=True)
  assert set(d.columns)==set(mapping.loc[mapping.Dataset.eq('ccRCC4')&mapping.TN.eq(tissue),'MetabID'])
 anno=pd.read_excel(x,sheet_name='metanno',index_col=0)
 anno.index=anno.index.astype(str).str.strip()
 oldpath=R/'results/tables/cohort_effects.tsv';old=pd.read_csv(oldpath,sep='\t');old=old[old.dataset.isin(['ccRCC3','ccRCC4'])]
 cancerpath=R/'results/tables/cancer_effects.tsv';frozen=pd.read_csv(cancerpath,sep='\t');frozen=frozen[frozen.cancer.eq('ccRCC')]
 save(old,pub/'historical_cohort_effects_preserved.tsv');save(frozen,pub/'historical_cancer_effects_preserved.tsv')
 old4=old[old.dataset.eq('ccRCC4')];assert old4.feature_name.is_unique
 inv=pd.DataFrame({'metabolite_name':mat.index})
 inv['feature_id']=inv.metabolite_name.map(lambda z:'ccRCC4_FEATURE:'+__import__('hashlib').sha256(z.encode()).hexdigest()[:16])
 inv=inv.merge(old4[['feature_name','metabolite_key','identity_confidence']],left_on='metabolite_name',right_on='feature_name',how='left',validate='one_to_one')
 for col in ['H_HMDB','H_KEGG','H_MetabolonID','H_name']:
  inv[col]=inv.metabolite_name.map(anno[col])
 inv['retained_both_author_filters']=inv.metabolite_name.isin(retained)
 inv['status']=np.where(inv.retained_both_author_filters,'DONE','NOT_EVALUABLE')
 inv['reason']=np.where(inv.retained_both_author_filters,'NA','not_retained_in_both_author_tissue_filters')
 save(inv,pub/'feature_inventory.tsv')
 summary=dict(status='DONE',ccRCC4_mapped_specimens=76,excluded_nonclear_or_discordant_specimens=len(excluded),excluded_patients=len(bad),
  ccRCC4_clear_cell_specimens=len(m),tumor_specimens=int(m.TN.eq('Tumor').sum()),normal_specimens=int(m.TN.eq('Normal').sum()),
  author_patients=m.SUBJECT_ID.nunique(),tumor_patients=m[m.TN.eq('Tumor')].SUBJECT_ID.nunique(),normal_patients=m[m.TN.eq('Normal')].SUBJECT_ID.nunique(),paired_patients=len(paired),
  original_matrix_specimens=len(mat.columns),unmapped_histology_or_rna_link_not_extended=len(mat.columns)-76,
  total_features=len(mat),retained_features=len(retained),historical_cancer_effects=len(frozen),historical_cancer_q005=int(frozen.effect_fdr.lt(.05).sum()),
  ccRCC3_status='DONE_SEPARATE_RUN',ccRCC3_reason='Resolved with MetaData_M4 in 20260925T151000Z_ccrcc3_v1;not pooled',
  treatment_scope='mixed treated and untreated, mostly immune-checkpoint treated;not a treatment-effect design',
  original_statistics_modified=False,observed_values_equal_between_sheets=True)
 spec=dict(version=V,code_commit=commit,cancer='ccRCC',cohort='ccRCC4',
  inclusion='MasterMapping Dataset=ccRCC4;both MasterMapping Histology and original-study MetaData_M5 HISTOLOGY equal Clear Cell;exclude whole patient on any non-clear-cell or discordant specimen',
  supersedes='20260925T143000Z_discovery_v1 ccRCC4 numeric results;initial CAMP-only histology mislabeled2patients;retain initial version as superseded history',
  patient_id='original OrigScale SUBJECT_ID, explicitly supplied by authors',
  analysis_unit='author patient; arithmetic mean of author processed values over all eligible regions within patient and tissue',
  discovery_scope='all features retained in both original author tissue-filtered sheets, before any significance selection; excluded feature inventory retained',
  available_sensitivity='patient retained for a feature only if ALL included regions in BOTH tissues have finite author data values; no changing regional mixture',
  transformation='none for metabolomics;author processed scale;no new imputation or batch adjustment',
  n_min=8,paired_patients=len(paired),bootstrap=4000,master_seed=20260925,
  test='signed rank;<=16 nonzero exact exhaustive signs;17-29 Monte Carlo99999;>=30 tie and continuity corrected normal;zero_method=wilcox',
  test_families=['METAB_PAIRED_PRIMARY','METAB_PAIRED_AVAILABLE'],BH='all evaluable features separately per family',
  workpool='primary P<0.05 exploratory, q separately retained',
  limitations=['small n','treated cohort','within-patient region mean estimand','same-cohort selection','cross-cohort independence not established','genotype identity not checked'],
  software=dict(python=platform.python_version(),pandas=pd.__version__,numpy=np.__version__),
  source_paper='https://doi.org/10.1038/s42255-023-00817-8')
 js(pub/'analysis_spec.json',spec);js(pub/'validation.json',summary)
 save(pd.DataFrame([dict(audit_item=k,value=v,status='NEEDS_REVIEW' if k.startswith('ccRCC3') else 'DONE') for k,v in summary.items()]),pub/'sample_audit_summary.tsv')
 treatment=m[['SUBJECT_ID','TREATMENT']].drop_duplicates().TREATMENT.value_counts().rename_axis('treatment').reset_index(name='n_patients');save(treatment,pub/'treatment_coverage.tsv')
 sources=[source,meta5path,mp,path,oldpath,cancerpath,Path(__file__),Path(__file__).with_name('ccrcc_stats_v1.py')]
 save(pd.DataFrame([dict(source_id=p.name,path_or_url=str(p),sha256=sha(p),access_scope='SERVER_ONLY_INPUT' if p.suffix=='.xlsx' or p==mp else 'SUMMARY_OR_CODE') for p in sources]),pub/'source_manifest.tsv')
 (pub/'README_CN.md').write_text(f'''# ccRCC样本设计与历史结果\n\n问题：能否按BRCA主线建立患者独立的ccRCC分析？\n\n输入：CAMP作者矩阵、MasterMapping和ccRCC4原始工作簿；逐样本表只留server165。\n\n实际结果：ccRCC4原76标本中排除2个髓质癌、1个嫌色癌、2个未分类癌标本对应的全部患者后，纳入{len(m)}个透明细胞癌标本、{m.SUBJECT_ID.nunique()}位作者患者、{len(paired)}对患者。多个区域先在患者×组织内等权取均值。全矩阵{len(mat)}特征，作者两组织均保留{len(retained)}项。原癌种冻结{len(frozen)}项及其P/q原样保存。\n\n新手解释：标本数不是患者数；同一患者的多个区域不能当作多位患者增加样本量。\n\n限制：治疗背景混合，主要为免疫治疗；不是治疗效果比较。ccRCC3患者连接已在另一运行解决，分别计算。本版本替代首版ccRCC4统计；首版保留但不进入主结论。矩阵内未被作者跨组学表确认的标本不擅自扩入。\n\n当前决定：按已冻结设计执行ccRCC4患者均值配对分析；保留ccRCC3状态。\n\n下一步：全量发现、直接关系、患者关联、RNA背景、全候选单细胞来源。\n\n复现：python ccrcc4_discovery_v2.py --out 新服务器运行目录 --code-commit 提交号；需预先放入source/ccRCC4_original.xlsx。\n''',encoding='utf-8')
 checksums(pub);print('AUDIT',json.dumps(summary),flush=True)
 p=out/'public/04_ROBUSTNESS';p.mkdir(exist_ok=False);js(p/'analysis_spec.json',spec)
 groups={t:m[m.TN.eq(t)].groupby('SUBJECT_ID').MetabID.apply(list).to_dict() for t in ['Tumor','Normal']}
 rows=[]
 for i,r in inv[inv.retained_both_author_filters].iterrows():
  name=r.metabolite_name
  means={t:np.array([mat.loc[name,groups[t][pid]].mean() for pid in paired]) for t in groups}
  avail=np.array([all(np.isfinite(raw.loc[name,groups[t][pid]]).all() for t in groups) for pid in paired])
  delta=means['Tumor']-means['Normal'];finite=np.isfinite(delta)
  for fam,mask in [('METAB_PAIRED_PRIMARY',finite),('METAB_PAIRED_AVAILABLE',finite&avail)]:
   d=dict.fromkeys(PREFIX,np.nan);d.update(cancer='ccRCC',cohort='ccRCC4',stage_id='04_ROBUSTNESS',run_id=out.name,analysis_version=V,analysis_type='patient_region_mean_paired',metabolite_key=r.metabolite_key,metabolite_name=name,gene='NA',unit='author_patient',test_family=fam,source_id='CAMP_ccRCC4_original_SUBJECT_ID',feature_id=r.feature_id,n_pairs_total=len(paired),n_both_available=int(avail.sum()),mask_semantics=spec['available_sensitivity'])
   d.update(paired_stats(delta[mask],fam,r.feature_id));d['effect_type']='mean_paired_difference_author_processed_scale';rows.append(d)
  if len(rows)%100==0:print('FEATURES',len(rows)//2,flush=True)
 d=adjust(pd.DataFrame(rows));assert not d.duplicated(['feature_id','test_family']).any();assert (d.n_up+d.n_down+d.n_equal).eq(d.n).all()
 for fam,fnm in [('METAB_PAIRED_PRIMARY','metabolite_paired.tsv'),('METAB_PAIRED_AVAILABLE','metabolite_available_sensitivity.tsv')]:save(d[d.test_family.eq(fam)],p/fnm)
 primary=d[d.test_family.eq('METAB_PAIRED_PRIMARY')].copy();primary['majority_direction_fraction']=primary[['up_fraction','down_fraction']].max(axis=1)
 primary['direction_discordant']=np.sign(primary.mean_delta).ne(np.sign(primary.n_up-primary.n_down))|np.sign(primary.mean_delta).ne(np.sign(primary.median_delta))
 primary=primary.sort_values(['majority_direction_fraction','p_value','feature_id'],ascending=[False,True,True])
 save(primary,p/'metabolite_direction_ranking.tsv');save(primary[primary.p_value.lt(.05)],p/'metabolite_workpool.tsv')
 val=dict(status='DONE',features=len(primary),paired_patients=len(paired),families=d.groupby('test_family').apply(lambda z:dict(planned=len(z),evaluable=int(z.p_value.notna().sum()),p005=int(z.p_value.lt(.05).sum()),q005=int(z.q_value.lt(.05).sum()))).to_dict(),BH_independently_verified=True,direction_counts_verified=True,original_statistics_modified=False)
 js(p/'validation.json',val);save(pd.read_csv(pub/'source_manifest.tsv',sep='\t'),p/'source_manifest.tsv')
 (p/'README_CN.md').write_text(f'''# ccRCC4患者配对代谢物发现\n\n问题：同患者肿瘤与正常组织的代谢物是否有变化？\n\n输入范围：经明确SUBJECT_ID、透明细胞病理筛选后的{len(paired)}对；{len(primary)}个作者两组织保留特征，多个区域先取患者内均值。\n\n实际结果：主分析P<0.05为{int(primary.p_value.lt(.05).sum())}项，q<0.05为{int(primary.q_value.lt(.05).sum())}项。详细检验族计数见validation.json。\n\n新手解释：P阈值用于探索性映射，q另列；患者内多个区域不增加独立样本量，效应不是浓度倍数。\n\n限制：小队列且多数有治疗背景；作者可用值敏感性要求患者全部纳入区域都有值，覆盖不足不算阴性。ccRCC3未合并。\n\n当前决定：按主分析P<0.05进入直接生化映射；完整结果及不显著背景保留。\n\n下一步：当前全关系患者关联、RNA背景和全基因细胞来源。\n\n复现：见code/ccrcc4_discovery_v2.py和analysis_spec.json。\n''',encoding='utf-8')
 checksums(p);(out/'.running').rename(out/'.discovery_done');print('DISCOVERY',json.dumps(val),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--code-commit',required=True);a=p.parse_args();main(a.out,a.code_commit)

