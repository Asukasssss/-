"""Server-only targeted paired-change analysis; no patient tables exported."""
import argparse, json, platform
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import prad_patient_v1 as engine

ap=argparse.ArgumentParser()
ap.add_argument('--out',type=Path,required=True)
ap.add_argument('--source-run',type=Path,required=True)
ap.add_argument('--code-commit',required=True)
a=ap.parse_args();out=a.out;src=a.source_run
assert out.is_dir()
with (out/'.running').open('x') as f:f.write('prad_slc6a6_paired_delta_v1')
pub=out/'public/03_PATIENT';pub.mkdir(parents=True,exist_ok=False)
private=out/'private';private.mkdir()
pp=src/'private/RNA_pairs_private.tsv';mp=src/'private/tumor_multiomics_map_private.tsv'
p=pd.read_csv(pp,sep='\t',dtype=str)
m=pd.read_csv(mp,sep='\t',dtype=str)
assert p['case'].is_unique and len(p)==43
for col in ['MetabID_tumor','MetabID_normal','RNAID_tumor','RNAID_normal']:assert p[col].is_unique
assert (p.Identifier_tumor==p.Identifier_normal).all()
rp=engine.SRC/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/m.RNAFile.iloc[0]
xp=engine.SRC/'processed_metabolomics/PreprocessedData_PRAD.xlsx'
rna=pd.read_csv(rp,index_col=0);met=pd.read_excel(xp,sheet_name='data_imputed',index_col=0);raw=pd.read_excel(xp,sheet_name='data',index_col=0)
for z in [met,raw]:z.index=z.index.astype(str).str.strip();assert z.index.is_unique
assert rna.index.is_unique
x=rna.loc['SLC6A6',p.RNAID_tumor].to_numpy(float)-rna.loc['SLC6A6',p.RNAID_normal].to_numpy(float)
y=met.loc['taurine',p.MetabID_tumor].to_numpy(float)-met.loc['taurine',p.MetabID_normal].to_numpy(float)
available=np.isfinite(raw.loc['taurine',p.MetabID_tumor].to_numpy(float))&np.isfinite(raw.loc['taurine',p.MetabID_normal].to_numpy(float))
assert np.isfinite(x).all() and np.isfinite(y).all()
assert ((x<0).sum(),(y<0).sum())==(32,36)
engine.save(pd.DataFrame(dict(case=p['case'],delta_RNA=x,delta_taurine=y,batch=p.Identifier_tumor,available=available,CAPT_concordant=p.CAPT_concordant)),private/'paired_changes_private.tsv')
spec=dict(version='prad_slc6a6_paired_delta_v1',code_parent_commit=a.code_commit,target='taurine|SLC6A6',selection='user-directed post-hoc single relationship; P<0.05 exploratory',unit='43 unique author case tumor-normal pairs',delta='tumor minus matched normal on unchanged author scales',RNA_scale='RMA BatchAdj continuous',metabolite_scale='author processed log2',families=['DELTA_PRIMARY','DELTA_AVAILABLE','DELTA_CAPT','DELTA_BATCH'],tests_per_family=1,q='single targeted test per family: q=p; not correction for earlier candidate selection',primary='Spearman correlation of paired deltas',permutations=9999,bootstrap=4000,seed_rule='reuse SHA256(prad_patient_v1|20260922|family|taurine|SLC6A6)',sensitivity='available original author data;36 CAPT-concordant pairs;batch rank residualization with within-batch permutations;leave-one-pair-out rho',direction='strict >0/<0;exact zero separately;no pseudocount or threshold',plots='unlabelled derived scatter;patient values and identifiers remain server-only',software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__))
engine.js(pub/'analysis_spec.json',spec)
rows=[];counts=[]
for family,keep,partial in [('DELTA_PRIMARY',np.ones(len(p),bool),False),('DELTA_AVAILABLE',available,False),('DELTA_CAPT',p.CAPT_concordant.eq('True').to_numpy(),False),('DELTA_BATCH',np.ones(len(p),bool),True)]:
 d=engine.base(out,family,'SLC6A6');d.update(metabolite_name='taurine',metabolite_key='KEGG:C00245',analysis_version=spec['version'])
 if family=='DELTA_AVAILABLE' and available.all():
  stats_result=dict(primary_stats)
  stats_result['reason']='identical_inputs_to_primary;reused_statistics'
 else:stats_result=engine.assoc(x[keep],y[keep],p.Identifier_tumor.to_numpy()[keep],family,'taurine|SLC6A6',partial)
 if family=='DELTA_PRIMARY':primary_stats=dict(stats_result)
 d.update(stats_result)
 d['n_pairs']=int(keep.sum())
 d['unit']='author_case_pair';d['ci_estimand']='paired_delta_rank_residual_correlation' if partial else 'paired_delta_Spearman_rho'
 d['p_method']=d['p_method'].replace('case','pair');d['ci_method']=d['ci_method'].replace('case','pair')
 d['q_value']=d['p_value'];d['family_n_evaluable']=1;rows.append(d)
 for sx in [-1,0,1]:
  for sy in [-1,0,1]:counts.append(dict(analysis_type=family,RNA_direction={-1:'down',0:'equal',1:'up'}[sx],taurine_direction={-1:'down',0:'equal',1:'up'}[sy],n=int(((np.sign(x[keep])==sx)&(np.sign(y[keep])==sy)).sum()),n_pairs=int(keep.sum())))
results=pd.DataFrame(rows);engine.save(results,pub/'results.tsv');engine.save(pd.DataFrame(counts),pub/'direction_counts.tsv')
primary=rows[0]
loo=np.array([stats.spearmanr(np.delete(x,i),np.delete(y,i)).statistic for i in range(len(x))])
engine.js(pub/'validation.json',dict(status='DONE',unique_case_pairs=43,CAPT_pairs=int(p.CAPT_concordant.eq('True').sum()),available_pairs=int(available.sum()),known_marginal_counts_reproduced=True,independent_Spearman_check=bool(np.isclose(stats.spearmanr(x,y).statistic,primary['rho'])),leave_one_out_rho_min=float(loo.min()),leave_one_out_rho_max=float(loo.max()),original_inputs_modified=False,independent_cohort=False))
fig,ax=plt.subplots(figsize=(7,6));ax.scatter(x,y,s=42,color='#287f9e',alpha=.8,edgecolors='white',linewidth=.5)
ax.axhline(0,color='#888888',lw=.8);ax.axvline(0,color='#888888',lw=.8)
ax.set_xlabel('SLC6A6 RNA change (tumor - matched normal)\nAuthor RMA BatchAdj scale')
ax.set_ylabel('Taurine change (tumor - matched normal)\nAuthor processed log2 scale')
ax.set_title(f"PRAD paired changes | n=43\nSpearman rho={primary['rho']:.3f}; permutation P={primary['p_value']:.4f}")
ax.spines[['top','right']].set_visible(False);fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(pub/f'paired_delta_scatter.{ext}',dpi=200)
plt.close(fig)
ct=pd.DataFrame(counts);tab=ct[ct.analysis_type.eq('DELTA_PRIMARY')]
report='\n'.join(f"- {r['analysis_type']}: n={r['n_pairs']}, rho={r['rho']:.6f}, P={r['p_value']:.6f}, bootstrap95%CI=[{r['ci_lower']:.6f},{r['ci_upper']:.6f}]" for r in rows)
signs='\n'.join(f'- RNA {r.RNA_direction}, taurine {r.taurine_direction}: {r.n}/43' for r in tab.itertuples())
(pub/'README_CN.md').write_text(f'''# PRAD 牛磺酸—SLC6A6 配对变化

本轮问题：同一患者 SLC6A6 的肿瘤-正常变化越低，牛磺酸变化是否也越低？

输入与范围：43对作者case配对；沿用作者处理尺度，不重新填补或变换；患者数据留服务器。单一关系为用户指定、同队列发现后的探索性补充。

实际结果：
{report}

方向组合：
{signs}

新手解释：每点为一个配对病例；横纵轴均为肿瘤减正常。正相关意味着两种变化量倾向同向排列。共同降低人数多，不自动证明下降幅度相关，尤其两者各自已有较高降低比例。零值单列，不剔除以提高一致性。

限制/反证：P<0.05为本轮展示标准，单关系各检验族q=P，不消除历史筛选偏倚。配对相关不证明因果、酶活或通量，也不能排除细胞组成、批次等影响。CAPT字段有7对不一致，36对一致样本单列敏感性。病例身份沿用作者case，非基因型验证。多个敏感性分析不是独立复现。原作者处理尺度下差值不是浓度倍数。

当前决定：全部结果保留，不按显著性选择图或样本。

下一步：独立样本与细胞组成控制尚未完成。

复现命令：python prad_slc6a6_paired_delta_v1.py --out NEW_RUN --source-run ORIGINAL_RUN --code-commit PARENT_COMMIT；先建立NEW_RUN并复制本脚本及prad_patient_v1.py。
''',encoding='utf-8')
engine.save(pd.DataFrame([dict(source_id=f.name,path_or_url=str(f),sha256=engine.sha(f)) for f in [pp,mp,rp,xp,Path(__file__),Path(engine.__file__)]]),pub/'source_manifest.tsv')
engine.save(pd.DataFrame([dict(file=f.name,sha256=engine.sha(f)) for f in sorted(pub.iterdir())]),pub/'checksums.tsv')
(out/'.running').rename(out/'COMPLETE')
print(report);print(signs);print('LOO',loo.min(),loo.max())
