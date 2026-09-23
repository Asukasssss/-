"""Paired donor-label inference on frozen single-cell aggregates. Server only."""
from pathlib import Path
import sys,json,hashlib,itertools,platform
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata,wilcoxon
from statsmodels.stats.multitest import multipletests
BASE=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A')
R=Path(sys.argv[1]); P=R/'public'; D=R/'private'
assert (R/'.running').read_text().strip()=='epithelial_paired20_v1'
P.mkdir();D.mkdir()
genes=sorted(pd.read_csv(R/'scope.tsv',sep='\t').gene.unique())
assert len(genes)==20
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
spec=dict(version='epithelial_paired20_v1',genes=genes,input_commit='ccc91bd',unit='publisher_donor_id',comparison='Malignant_epithelial minus Nonmalignant_epithelial within same donor',reference='author-annotated nonmalignant cells within tumor;not healthy breast',min_cells_each_type=20,min_paired_donors=6,min_max_group_detection=.01,test='two-sided exact sign enumeration of Wilcoxon signed ranks;zero differences excluded;average tied ranks',effect='mean paired difference in mean cellular log1p(count/library*10000);not log fold change',ci='10000 paired-donor bootstrap percentile pointwise95;not simultaneous',families='Wu ALL20 main;Wu author Naive20 sensitivity;Pal20 not evaluable. BH20 with unavailable P=1 internally;Holm20 also reported. Display P stays NA for unavailable',seed=20260923,post_selection_exploratory=True,cell_level_tests=False,old_statistics_modified=False,unpaired_test=False,annotation='frozen author celltype_major;no new CNV calls',independence_limit='source donor labels;not genotype checked',software={'python':platform.python_version(),'scipy':scipy.__version__,'numpy':np.__version__,'pandas':pd.__version__})
(P/'analysis_spec.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
inputs=[]; frames={}
for cohort in ['Wu2021','Pal2021_reprocessed']:
    parts=[]
    for run in ['20260919T140105Z_scRNA117_v1','20260922T023000Z_sc39_v1']:
        p=BASE/run/'private'/(cohort+'_donor_profiles.tsv')
        inputs.append(dict(path=str(p),sha256=sha(p),kind='frozen_private_donor_profile'))
        x=pd.read_csv(p,sep='\t');parts.append(x[x.gene.isin(genes)])
    x=pd.concat(parts,ignore_index=True)
    assert set(x.gene)==set(genes)
    assert not x.duplicated(['donor','celltype','gene']).any(),'Multiple strata must be resolved before inference'
    frames[cohort]=x
save(pd.DataFrame(inputs),P/'source_manifest.tsv')
rows=[];cover=[];private=[];verification=[]
M='Malignant_epithelial';N='Nonmalignant_epithelial'
for cohort,x in frames.items():
  for partition in (['ALL','AUTHOR_NAIVE'] if cohort=='Wu2021' else ['ALL']):
    d=x.copy()
    if partition=='AUTHOR_NAIVE': d=d[d.treatment.eq('Naïve')]
    for gene in genes:
      z=d[d.gene.eq(gene)&d.celltype.isin([M,N])&d.n_cells.ge(20)]
      pivot=z.pivot(index='donor',columns='celltype',values='mean_log1p10k').reindex(columns=[M,N]); paired=pivot.dropna()
      n=len(paired)
      det=z.pivot(index='donor',columns='celltype',values='detection_fraction').reindex(index=paired.index,columns=[M,N])
      cov=dict(cohort=cohort,partition=partition,gene=gene,n_malignant=int(pivot[M].notna().sum()),n_nonmalignant=int(pivot[N].notna().sum()),n_paired=n)
      cover.append(cov)
      row=dict(cancer='BRCA',cohort=cohort,stage_id='06_EXTERNAL',run_id=R.name,analysis_version=spec['version'],analysis_type=partition,metabolite_key='NA',metabolite_name='NA',gene=gene,unit='paired_publisher_donor_label',n=n,n_reference=n,effect_type='mean_paired_mean_log1p10k_difference',effect=np.nan,ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,q_value=np.nan,test_family=cohort+'_'+partition+'_planned20',family_n_evaluable=0,status='NOT_EVALUABLE',reason='',source_id='GSE176078' if cohort=='Wu2021' else 'Pal2021_reprocessed',holm20=np.nan,positive_pairs=0,negative_pairs=0,equal_pairs=0,mean_malignant=np.nan,mean_nonmalignant=np.nan,mean_detection_malignant=np.nan,mean_detection_nonmalignant=np.nan,median_difference=np.nan)
      if n:
        delta=paired[M].values-paired[N].values
        row.update(effect=float(delta.mean()),mean_malignant=float(paired[M].mean()),mean_nonmalignant=float(paired[N].mean()),mean_detection_malignant=float(det[M].mean()),mean_detection_nonmalignant=float(det[N].mean()),positive_pairs=int((delta>0).sum()),negative_pairs=int((delta<0).sum()),equal_pairs=int((delta==0).sum()),median_difference=float(np.median(delta)))
        for donor in paired.index:
            private.append(dict(cohort=cohort,partition=partition,gene=gene,donor=donor,malignant=paired.loc[donor,M],nonmalignant=paired.loc[donor,N],difference=paired.loc[donor,M]-paired.loc[donor,N]))
      if n<6:row['reason']='fewer_than_6_paired_donors;no_unpaired_fallback'
      elif max(det[M].mean(),det[N].mean())<.01:row['reason']='both_groups_mean_detection_below_1percent'
      else:
        nz=delta[delta!=0]; ranks=rankdata(abs(nz));obs=abs(np.dot(np.sign(nz),ranks))
        signs=np.array(list(itertools.product([-1,1],repeat=len(nz))))
        p=float(np.mean(abs(signs@ranks)>=obs-1e-12)) if len(nz) else 1.
        rng=np.random.default_rng(20260923+genes.index(gene))
        b=delta[rng.integers(0,n,size=(10000,n))].mean(axis=1);ci=np.quantile(b,[.025,.975])
        row.update(status='DONE',reason='',p_value=p,ci_lower=float(ci[0]),ci_upper=float(ci[1]))
        if len(nz)==n and len(np.unique(abs(nz)))==n:
            p2=float(wilcoxon(delta,alternative='two-sided',method='exact').pvalue)
            assert abs(p-p2)<1e-12
            verification.append(dict(gene=gene,partition=partition,enumerated_p=p,scipy_p=p2))
      rows.append(row)
out=pd.DataFrame(rows)
for family,idx in out.groupby('test_family').groups.items():
    good=out.loc[idx,'status'].eq('DONE');ps=out.loc[idx,'p_value'].fillna(1).values
    for col,method in [('q_value','fdr_bh'),('holm20','holm')]:
        adj=multipletests(ps,method=method)[1];out.loc[np.array(list(idx))[good.values],col]=adj[good.values]
    out.loc[idx,'family_n_evaluable']=int(good.sum())
out['conclusion']=np.where(out.status.ne('DONE'),'NOT_EVALUABLE',np.where(out.q_value.lt(.05)&out.effect.gt(0),'HIGHER_FDR',np.where(out.q_value.lt(.05)&out.effect.lt(0),'LOWER_FDR','NOT_FDR_SUPPORTED')))
save(out,P/'results.tsv');save(pd.DataFrame(cover),P/'coverage.tsv');save(pd.DataFrame(private),D/'paired_measurements.tsv');save(pd.DataFrame(verification),P/'independent_p_check.tsv')
counts=out.groupby(['cohort','analysis_type','conclusion']).size().reset_index(name='count');save(counts,P/'summary_counts.tsv')
valid=dict(status='PASS',scope20=True,rows=len(out),duplicate_donor_gene_type=False,exact_p_scipy_checks=len(verification),min_cells=20,private_measurements_exported=False,source_profile_checksums=True,patient_identity_reverified=False,source_matrix_reread=False,visual_review='PENDING')
(P/'validation.json').write_text(json.dumps(valid,indent=2),encoding='utf-8')
(P/'README_CN.md').write_text('''# 专题20基因：恶性与非恶性上皮的供者配对比较

## 本轮问题
检验恶性上皮是否高于同供者非恶性上皮；不预设升高。
## 输入与范围
复用Wu/Pal原117及新增39供者汇总，限定专题20基因。每供者每类至少20细胞，至少6位配对供者。Wu实际8位；作者Naive标签敏感性7位。Pal当前子集无非恶性上皮，20项不可评估，不构成阴性或复现。
## 实际结果
见results.tsv完整主分析、作者Naive敏感性和Pal缺项。表内效应来自两类均有覆盖的相同供者，不能与原点图不同供者组成的平均排序混用。
## 新手解释
每位供者产生一个恶性减非恶性表达差，再进行双侧精确符号秩检验。效应为归一化单细胞log1p表达的供者均值差，不是浓度倍数或log2FC。方向人数不是独立显著人数。
## 限制与反证
非恶性上皮来自肿瘤标本内作者注释，不是健康人或癌旁正常组织。沿用作者供者ID，未基因型验证。未调整上皮亚群组成、亚型或治疗；作者Naive仅作标签敏感性，不新解释治疗时间。8位覆盖有限，且本专题已查看表达图，为后选择探索；新BH20/Holm20不覆盖整个项目历史选择。每个区间是点区间，不是全20项同时区间。低检出不证明无功能。
## 当前决定
按供者检验结果更新表达证据；不能单凭表达确定靶点或正常组织安全性。
## 下一步
本轮到差异检验与汇总停止；无新聚类、无细胞级伪重复、无外部资源下载。
## 复现
在server165新运行目录创建独占.running，内容epithelial_paired20_v1；复制scope.tsv及脚本，python3 brca_epithelial_paired20_v1.py RUN_DIR。完整供者数值仅存private。来源：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078 。
''',encoding='utf-8')
print(counts.to_string(index=False));print(out[out.analysis_type.eq('ALL')&out.cohort.eq('Wu2021')][['gene','n','effect','positive_pairs','p_value','q_value','holm20','conclusion']].to_string(index=False))
(R/'.running').rename(R/'.completed')
