"""Descriptive pseudobulk CPM fold; not exponentiation of prior log-mean contrast."""
from pathlib import Path
import sys,json,hashlib
import numpy as np,pandas as pd
B=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A')
R=Path(sys.argv[1]);assert (R/'.running').read_text().strip()=='fold20_v1'
P=R/'public';P.mkdir();D=R/'private';D.mkdir()
previous=B/'20260923T040012Z_epithelial_paired20_v1'
prior=pd.read_csv(previous/'public/results.tsv',sep='\t');prior=prior[prior.cohort.eq('Wu2021')&prior.analysis_type.eq('ALL')]
paired=pd.read_csv(previous/'private/paired_measurements.tsv',sep='\t');paired=paired[paired.cohort.eq('Wu2021')&paired.partition.eq('ALL')]
paths=[B/n/'private/Wu2021_donor_profiles.tsv' for n in ['20260919T140105Z_scRNA117_v1','20260922T023000Z_sc39_v1']]
x=pd.concat([pd.read_csv(p,sep='\t') for p in paths]);M='Malignant_epithelial';N='Nonmalignant_epithelial'
rows=[];priv=[]
for r in prior.itertuples():
    donors=paired.loc[paired.gene.eq(r.gene),'donor'];assert len(donors)==8
    z=x[x.gene.eq(r.gene)&x.donor.isin(donors)&x.celltype.isin([M,N])].copy()
    assert len(z)==16 and not z.duplicated(['donor','celltype']).any() and z.n_cells.ge(20).all()
    z['CPM']=z.raw_count/z.all_gene_library_sum*1e6
    t=z.pivot(index='donor',columns='celltype',values='CPM');assert len(t)==8 and not t.isna().any().any()
    a,b=t[M].mean(),t[N].mean();rat=t[M].div(t[N].replace(0,np.nan))
    row=dict(gene=r.gene,n=8,mean_malignant_pseudobulk_CPM=a,mean_nonmalignant_pseudobulk_CPM=b,ratio_of_equal_donor_means=a/b if b>0 else np.nan,median_within_donor_ratio=rat.median(),n_ratio_evaluable=int(rat.notna().sum()),n_zero_reference=int(t[N].eq(0).sum()),CPM_higher_pairs=int(t[M].gt(t[N]).sum()),prior_logmean_higher_pairs=r.positive_pairs,prior_logmean_p=r.p_value,prior_logmean_q=r.q_value,prior_status=r.status,status='DONE' if b>0 else 'NOT_EVALUABLE',reason='' if b>0 else 'zero_mean_reference_no_pseudocount')
    rows.append(row)
    for donor in t.index:priv.append(dict(gene=r.gene,donor=donor,malignant_CPM=t.loc[donor,M],nonmalignant_CPM=t.loc[donor,N],ratio=rat.loc[donor]))
out=pd.DataFrame(rows);out.to_csv(P/'fold20.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(priv).to_csv(D/'paired_CPM.tsv',sep='\t',index=False,na_rep='NA')
spec=dict(method='per donor/type pseudobulk sum gene counts / sum all-gene library counts * 1e6;ratio of equal-donor mean CPM',n=8,genes=prior.gene.tolist(),new_P_q=0,no_pseudocount=True,prior_p_q='separate historical columns;test mean cellular log1p10k,not this CPM ratio',limits='relative RNA abundance,not absolute RNA molecules,protein,activity or malignant selectivity;cell library-size weights differ from prior logmean estimand',median_ratio='secondary;zero reference excluded explicitly')
(P/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
paths += [previous/'private/paired_measurements.tsv',previous/'public/results.tsv',Path(__file__)]
pd.DataFrame([dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths]).to_csv(P/'source_manifest.tsv',sep='\t',index=False)
(P/'README_CN.md').write_text('''# 8位供者上皮表达倍数补充
## 本轮问题
回答恶性上皮相对非恶性上皮表达多少倍。
## 输入与范围
与上一轮相同8位供者、20基因；每类至少20细胞，沿用作者注释。完整供者数值仅存服务器。
## 实际结果
见fold20.tsv。主倍数为供者等权的两组平均CPM之比，不是逐供者倍数的平均；另保留可计算的逐供者倍数中位数及参照为0人数。
## 新手解释
每位供者每类细胞汇总基因计数，除以该类全部基因总计数得到CPM，然后对8位供者等权平均；恶性均值除以非恶性均值。3倍表示约为参照的3倍，即高约200%。
## 限制与反证
这是相对RNA丰度描述，不是绝对分子数、蛋白或酶活。不是将原log均值差指数化，归一化和细胞权重与原检验不同。旧P/q只在明确前缀列保留，不作为这个CPM倍数的检验。低检出基因即使可算倍数也不因此升级。
## 当前决定与下一步
补充效应大小；不改变上一轮0项BH20支持的结论。此处停止。
## 复现
server165独占目录.running内容fold20_v1，python3 brca_epithelial_fold20_v1.py RUN_DIR。
''',encoding='utf-8')
(P/'validation.json').write_text(json.dumps(dict(status='PASS',same8_donors_per_gene=True,no_duplicates=True,min_cells20=True,new_P_q=0,private_export=False),indent=2))
(P/'.gitattributes').write_bytes(b'* -text\n')
pd.DataFrame([dict(path=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(P.iterdir())]).to_csv(P/'checksums.tsv',sep='\t',index=False)
print(out[out.gene.isin(['LYPLA1','LYPLA2','ABHD12','LPCAT1','PCYT2'])].to_string(index=False))
(R/'.running').rename(R/'.completed')
