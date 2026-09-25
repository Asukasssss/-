"""Donor-paired SLC6A6 cancer versus adjacent benign cell-class comparison."""
import argparse, hashlib, json, platform
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--source-run',type=Path,required=True);ap.add_argument('--code-commit',required=True);a=ap.parse_args()
out=a.out
with (out/'.running').open('x') as f:f.write('prad_slc6a6_sc_paired_v1')
pub=out/'public/06_EXTERNAL';pub.mkdir(parents=True,exist_ok=False)
private=out/'private';private.mkdir()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA')
def js(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
src=a.source_run/'private/sc_donor_profiles_private.tsv'
allp=pd.read_csv(src,sep='\t');p=allp[allp.gene.eq('SLC6A6')].copy()
assert not p.duplicated(['partition','donor','celltype']).any()
assert set(p.partition)=={'cancer','adj_benign'}
spec=dict(version='prad_slc6a6_sc_paired_v1',parent_commit=a.code_commit,gene='SLC6A6',unit='matched depositor donor_id within cell class',source='same PRAD24 CELLxGENE dataset as source run',expression='reuse donor mean log1p(CP10K) from full raw-gene library normalization',minimum_cells_each_donor_tissue_celltype=20,minimum_pairs_for_inference=8,primary_family='all author same-label classes except Unassigned; epithelial subclasses retained',secondary_family='tumor malignant epithelium versus adjacent normal epithelium;different states explicitly labelled',primary_endpoint='mean per-cell normalized expression within donor and class',test='two-sided signed-rank exact sign-flip of nonzero differences;ties average ranks',bootstrap=4000,seed=20260925,display_threshold='P<0.05 exploratory;BH retained separately per family',detection_endpoint='descriptive only;no detection P tests',limitations=['same study as source ranking, not independent confirmation','no cell-level pseudoreplication','donor IDs and annotations inherited from depositor','mean normalized expression rather than count pseudobulk model','composition within broad cell classes and sample processing may remain confounded'])
js(pub/'analysis_spec.json',spec)
classes=sorted(set(p.celltype)-{'Unassigned'})
contrasts=[(c,c,'SAME_CLASS',c) for c in classes]+[('Epithelial_malignant','Epithelial_normal','EPITHELIAL_STATE','Malignant_vs_adjacent_normal_epithelium')]
rows=[];private_rows=[]
for tc,nc,family,label in contrasts:
 t=p[p.partition.eq('cancer')&p.celltype.eq(tc)&p.n_cells.ge(20)]
 n=p[p.partition.eq('adj_benign')&p.celltype.eq(nc)&p.n_cells.ge(20)]
 pair=t.merge(n,on='donor',suffixes=('_tumor','_normal'),validate='one_to_one')
 x=pair.mean_log1p_cp10k_tumor.to_numpy();y=pair.mean_log1p_cp10k_normal.to_numpy();delta=x-y;N=len(delta)
 pair['contrast']=label;private_rows.append(pair)
 r=dict(cancer='PRAD',cohort='PRAD24_CELLxGENE',stage_id='06_EXTERNAL',run_id=out.name,analysis_version=spec['version'],analysis_type='paired_celltype_expression',metabolite_key='NA',metabolite_name='NA',gene='SLC6A6',unit='paired_donor',n=N,n_reference=15,effect_type='mean_delta_donor_mean_log1p_CP10K',effect=float(delta.mean()) if N else np.nan,ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,q_value=np.nan,test_family=family,family_n_evaluable=0,status='NOT_EVALUABLE',reason='fewer_than8_eligible_donor_pairs',source_id='PRAD24_CELLxGENE',contrast=label,tumor_celltype=tc,adjacent_celltype=nc,n_up=int((delta>0).sum()),n_down=int((delta<0).sum()),n_equal=int((delta==0).sum()),mean_tumor=float(x.mean()) if N else np.nan,mean_adjacent=float(y.mean()) if N else np.nan,mean_detection_tumor=float(pair.detection_fraction_tumor.mean()) if N else np.nan,mean_detection_adjacent=float(pair.detection_fraction_normal.mean()) if N else np.nan)
 if N>=8:
  nz=delta[delta!=0];ranks=stats.rankdata(abs(nz));obs=float(np.sum(np.sign(nz)*ranks))
  assert len(nz)<=20
  bits=(np.arange(2**len(nz))[:,None]>>np.arange(len(nz)))&1
  null=(2*bits-1)@ranks
  pv=float(np.mean(abs(null)>=abs(obs)-1e-12))
  # Independent scipy exhaustive permutation check of signed-rank statistic.
  check=stats.wilcoxon(delta,zero_method='wilcox',method=stats.PermutationMethod(n_resamples=np.inf)).pvalue if len(nz) else 1.
  assert np.isclose(pv,check,atol=1e-12),(label,pv,check)
  seed=20260925+int(hashlib.sha256(label.encode()).hexdigest()[:6],16)
  rng=np.random.default_rng(seed);bs=delta[rng.integers(N,size=(4000,N))].mean(1)
  r.update(p_value=pv,ci_lower=float(np.quantile(bs,.025)),ci_upper=float(np.quantile(bs,.975)),status='DONE',reason='exploratory_same_study',seed=seed,n_nonzero=len(nz),exact_sign_patterns=len(null))
 rows.append(r)
d=pd.DataFrame(rows)
for family,ix in d.groupby('test_family').groups.items():
 ok=d.index.isin(ix)&d.p_value.notna();d.loc[ix,'family_n_evaluable']=int(ok.sum())
 if ok.any():d.loc[ok,'q_value']=multipletests(d.loc[ok,'p_value'],method='fdr_bh')[1]
save(d,pub/'results.tsv');save(pd.concat(private_rows,ignore_index=True),private/'paired_donor_values_private.tsv')
eligible=d[d.status.eq('DONE')].copy()
fig,ax=plt.subplots(figsize=(9,max(4,len(eligible)*.48)))
for i,(_,r) in enumerate(eligible.iterrows()):
 ax.plot([r.ci_lower,r.ci_upper],[i,i],color='#777777');ax.scatter(r.effect,i,color='#bb453f' if r.p_value<.05 else '#287f9e',zorder=3)
ax.set_yticks(range(len(eligible)));ax.set_yticklabels([f'{r.contrast} (n={r.n}, P={r.p_value:.3g})' for r in eligible.itertuples()],fontsize=9)
ax.axvline(0,color='#aaaaaa',lw=.8);ax.set_xlabel('Tumor - adjacent: donor mean log1p(CP10K)\nPoints: mean difference; lines: 95% paired bootstrap CI')
ax.set_title('SLC6A6 | paired donors within cell classes');ax.spines[['top','right']].set_visible(False);fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(pub/f'SLC6A6_paired_celltypes.{ext}',dpi=180)
plt.close(fig)
js(pub/'validation.json',dict(status='DONE',unique_donor_class_keys=True,tests_checked_against_scipy_exhaustive=True,total_contrasts=len(d),evaluable_contrasts=len(eligible),source_hash=sha(src),source_unchanged=True,patient_tables_exported=False))
lines='\n'.join(f'- {r.contrast}: n={r.n}, effect={r.effect:.6g}, P={r.p_value:.6g}, up/down/equal={r.n_up}/{r.n_down}/{r.n_equal}, status={r.status}' for r in d.itertuples())
(pub/'README_CN.md').write_text(f'''# SLC6A6 单细胞供者配对比较

本轮问题：同类细胞中，SLC6A6在肿瘤组织与癌旁组织是否不同？

输入与范围：沿用作者细胞标签及donor_id，复用前轮供者表达汇总；每位供者两种组织对应类别均至少20细胞，至少8对才做推断。所有同名细胞类别预先纳入；另列恶性上皮与癌旁正常上皮的状态比较，不能当作完全相同细胞状态。匹配仅按明确donor_id，不按顺序推断。

实际结果：
{lines}

新手解释：以供者为统计单位，每位供者先计算同类细胞平均表达，再比较肿瘤减癌旁。P来自双侧精确符号置换的Wilcoxon有符号秩检验；展示按P<0.05，q保留不作为入口。图中区间为平均变化的配对bootstrap区间，与秩检验不是同一统计量。

限制/反证：同一研究的补充分析，不是第二队列验证；低覆盖类别不可评估，不作阴性。表达是平均log1p(CP10K)，不是浓度倍数或转运活性。宽细胞类型内部的亚群构成变化、技术因素仍可影响结果。继承原CELLxGENE固定版本和其GEO链接不一致警告。

当前决定：按真实结果报告，不将组织RNA差异套用到某一细胞类型。

下一步：独立研究和亚群构成检查尚未完成。

复现命令：python prad_slc6a6_sc_paired_v1.py --out NEW_RUN --source-run ORIGINAL_RUN --code-commit PARENT；先建立新运行目录。
''',encoding='utf-8')
save(pd.DataFrame([dict(source_id=f.name,path_or_url=str(f),sha256=sha(f)) for f in [src,Path(__file__),a.source_run/'public/06_EXTERNAL/analysis_spec.json']]),pub/'source_manifest.tsv')
save(pd.DataFrame([dict(file=f.name,sha256=sha(f)) for f in sorted(pub.iterdir())]),pub/'checksums.tsv')
(out/'.running').rename(out/'COMPLETE')
print(d[['contrast','n','mean_tumor','mean_adjacent','effect','p_value','q_value','n_up','n_down','status']].to_string(index=False))
