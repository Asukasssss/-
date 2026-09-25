"""Nominal-P only report and figures; no FDR calculation or selection."""
import argparse,json,sys,hashlib,textwrap
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'runtime/report_dependencies'))
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();out=a.out
d=pd.read_csv(out/'results.tsv',sep='\t');o=pd.read_csv(out/'reactome_ORA.tsv',sep='\t');g=pd.read_csv(out/'reactome_GSEA.tsv',sep='\t');v=json.loads((out/'validation.json').read_text())
sel=(d.p_value<.05)&(d.log2FC.abs()>=.5)
assert sel.sum()==v['selected_high']+v['selected_low'] and (d.p_value<.05).sum()==v['gene_P_below_05']
assert np.allclose(hypergeom.sf(o.overlap-1,o.background_genes,o.pathway_tested_genes,o.selected_genes),o.p_value,rtol=1e-10,atol=1e-15)
assert d.q_value.isna().all() and 'padj' not in g and 'q_value' not in o and 'LYPLA1' not in set(d.gene)
# Compact repeated descriptive metadata for the repository's 5 MB file limit.
# Numerical values, gene universe and complete shared prefix remain unchanged.
d['analysis_type']='LYPLA1_high_low'
d['reason']='NOMINAL_P_ONLY'
d['test_family']='11058_gene_nominal_P'
for name,frame in [('results.tsv',d),('selected_genes.tsv',d[sel]),('all_genes_P05.tsv',d[d.p_value<.05])]:
    frame.to_csv(out/name,sep='\t',index=False,na_rep='NA',lineterminator='\n')
hits=[];lookup=d.set_index('gene')
for _,r in o[o.p_value<.05].iterrows():
    for gene in str(r.overlap_genes).split(';'):
        x=lookup.loc[gene];hits.append(dict(pathway_id=r.pathway_id,pathway=r.pathway,direction=r.direction,pathway_P=r.p_value,gene=gene,log2FC=x.log2FC,gene_P=x.p_value,depth_log2FC=x.depth_log2FC,depth_P=x.depth_p_value,donors_high_greater=x.donors_high_greater,donors_high_lower=x.donors_high_lower))
pd.DataFrame(hits).to_csv(out/'pathway_gene_links.tsv',sep='\t',index=False,lineterminator='\n')
plt.rcParams['font.family']='Microsoft YaHei'
fig,ax=plt.subplots(figsize=(9,6));colors=np.where(sel,np.where(d.log2FC>0,'#c86639','#397f9b'),'#bbc0c5')
ax.scatter(d.log2FC,-np.log10(d.p_value.clip(lower=1e-300)),c=colors,s=9,alpha=.65,rasterized=True)
ax.axhline(-np.log10(.05),c='gray',ls='--',lw=.8)
for x in [-.5,.5]:ax.axvline(x,c='gray',ls='--',lw=.8)
for _,r in d[sel].nsmallest(8,'p_value').iterrows():ax.annotate(r.gene,(r.log2FC,-np.log10(r.p_value)),fontsize=8,xytext=(4,4),textcoords='offset points')
ax.set(xlabel='log2倍数变化（LYPLA1高 / 低）',ylabel='−log10(原始P)',title=f'LYPLA1阳性CNA上皮｜供者内高低表达比较\nP<0.05且|log2FC|≥0.5：高组较高{v["selected_high"]}个，较低{v["selected_low"]}个')
fig.tight_layout();fig.savefig(out/'differential_P.png',dpi=180);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(17,7))
for ax,direction,title,color in zip(axes,['higher_in_LYPLA1_high','lower_in_LYPLA1_high'],['高组较高基因的通路富集','高组较低基因的通路富集'],['#c86639','#397f9b']):
    z=o[(o.direction==direction)&(o.p_value<.05)].nsmallest(8,'p_value').iloc[::-1]
    ax.barh(range(len(z)),-np.log10(z.p_value),color=color)
    ax.set_yticks(range(len(z)),['\n'.join(textwrap.wrap(t,32)) for t in z.pathway],fontsize=9)
    ax.set_title(title);ax.set_xlabel('−log10(原始P)')
    for i,(_,r) in enumerate(z.iterrows()):ax.text(-np.log10(r.p_value)+.025,i,f'n={r.overlap}',va='center',fontsize=8)
    ax.set_xlim(0,max(-np.log10(z.p_value))*1.18)
fig.suptitle('Reactome差异基因富集｜每个方向展示P最小的8项\nn为命中基因数；仅按原始P筛选，通路间可共享基因',fontsize=15)
fig.tight_layout(rect=[0,0,1,.91]);fig.savefig(out/'reactome_ORA_P.png',dpi=180);plt.close(fig)
top=d[sel].head(12);table='\n'.join(f'|{r.gene}|{r.log2FC:.3f}|{r.p_value:.3g}|{r.depth_p_value:.3g}|' for _,r in top.iterrows())
ot=o[o.p_value<.05].head(15);otab='\n'.join(f'|{r.pathway}|{"高组较高" if r.direction=="higher_in_LYPLA1_high" else "高组较低"}|{r.p_value:.3g}|{r.overlap_genes}|' for _,r in ot.iterrows())
txt=f'''# COAD：LYPLA1高低表达分组差异与通路（原始P口径）

## 本轮问题

按用户要求，不用FDR筛选或展示，以原始P<0.05探索LYPLA1高、低表达CNA上皮的伴随基因和通路。患者分组、原始P及效应均复用；仅改变筛选口径并为新基因列表重新做ORA。不把富集视为因果调控。

## 输入与范围

Uhlitz GSE166555原作者CNA肿瘤上皮；沿用来源冲突供者排除规则。4,477个CNA细胞中2,549个检出LYPLA1；仅在阳性细胞内，按各供者的log1p(CP10K)中位数分为高（严格大于中位数）、低（其余）组。每组至少20细胞，1位供者不足，最终8位供者、2,528细胞，高1,262、低1,266。1,928个零计数不混入低组。

全基因源计数21,854项，过滤后11,058项进入edgeR供者配对伪合并比较（TMM、稳健QL，`~patient+group`），LYPLA1自身排除。主筛选P<0.05，另保留|log2FC|≥0.5作为效应门槛；完整P<0.05清单也单独提供。深度敏感性加入每个供者组别的平均UMI对数，不作为挑选主模型的依据。方法见[Bioconductor](https://www.bioconductor.org/books/3.19/OSCA.multisample/multi-sample-comparisons.html)。

## 实际结果

{v['gene_P_below_05']}个基因P<0.05；其中{v['selected_high']+v['selected_low']}个同时达到效应阈值，高组较高{v['selected_high']}个、较低{v['selected_low']}个。深度敏感性模型中有{v['depth_selected']}个达到相同P及效应门槛，与主清单交集{v['primary_and_depth_selected']}个。主清单中{v['primary_selected_same_direction_depth_P05']}个在深度模型仍P<0.05且方向一致；这是同一数据的敏感性，不是复现。

![差异基因](differential_P.png)

|基因|高/低log2FC|原始P|深度模型P|
|---|---:|---:|---:|
{table}

[Reactome人通路](https://reactome.org/download-data)按实际可检验基因取交集后，15–500基因的通路共{v['pathways_tested']}条。ORA以上述131个基因为输入、高低方向分别分析，全部11,058个可检验基因为背景；两方向合计{v['ORA_P05']}个条目P<0.05。GSEA复用全基因排序，{v['GSEA_P05']}条P<0.05，高组端{v['GSEA_high_P05']}条、低组端{v['GSEA_low_P05']}条。两种方法不是独立证据，不相加计数。

![通路](reactome_ORA_P.png)

|ORA通路|输入方向|原始P|命中基因|
|---|---|---:|---|
{otab}

## 新手解释

ORA回答选出的差异基因是否较多落在某通路；GSEA使用全部基因的带方向排序，正NES偏高组、负NES偏低组。通路命中和leadingEdge基因并非同一概念，后者不要求每个基因单独P<0.05。

完整文件：`all_genes_P05.tsv`（全部945项）、`selected_genes.tsv`（加效应阈值131项）、`results.tsv`（11,058项含敏感性）、`reactome_ORA.tsv`、`reactome_GSEA.tsv`、`pathway_gene_links.tsv`（每条P<0.05 ORA通路的具体基因与各自统计）。统一结果表的q列为NA，仅维持跨癌字段，旧q不用于本版；历史原文件未覆盖。

## 限制/反证

本版全部按未校正P探索。很多条目只由2–4个共享基因支持，例如NMU/KISS1/PPBP会出现在多个GPCR条目，不能把这些条目当作不同机制的重复支持。GSEA低组端的多个翻译/核糖体条目共享大量RPL/RPS基因；不能据命名推断病毒感染、氨基酸缺乏或具体通路活性。

高组总UMI中位数9,502，低组13,027.5，分组与深度有联系；敏感性也不能排除全部混杂。部分基因仅少数供者检出，完整表保留供者差值方向数，不能只看最小P。肿瘤上皮内仍可能有克隆、周期、分化及背景RNA差异。CNA是作者RNA推断，不是逐细胞DNA认证；本分析是表达关联，未做LYPLA1干预或代谢测量。

## 当前决定

交付原始P口径的差异基因、通路及命中基因清单。保留全量与深度敏感性，不把131个基因或35个通路自动升级为LYPLA1下游靶点。

## 下一步

本批到差异与富集探索收口，先结合具体命中基因阅读；不自动扩展细胞通讯或拟时序。

## 复现命令

原始配对差异脚本`analyze_lypla1_highlow_v1.R`与分组参数已版本化；本次在server165新独占目录运行`review_lypla1_highlow_nominal_v2.py --out <新目录>`，复用旧效应、P和GSEA，按新基因列表计算ORA。本地`report_lypla1_highlow_nominal_v2.py --out <结果目录>`独立核查超几何P和计数后绘图。患者/细胞数据与Reactome源文件留服务器；只发布汇总。
'''
(out/'README_CN.md').write_text(txt,encoding='utf-8',newline='\n')
(out/'arithmetic_validation.json').write_text(json.dumps(dict(status='PASS',ORA_hypergeometric_P=True,selection_counts=True,no_new_FDR=True),indent=2)+'\n',encoding='utf-8')
files=['prepare_lypla1_highlow_v1.py','analyze_lypla1_highlow_v1.R','review_lypla1_highlow_nominal_v2.py','report_lypla1_highlow_nominal_v2.py','lypla1_highlow_v1.json']
(out/'code_manifest.json').write_text(json.dumps([dict(path='code/coad/'+f,sha256=hashlib.sha256((ROOT/'code/coad'/f).read_bytes()).hexdigest()) for f in files],indent=2)+'\n',encoding='utf-8')
for p in out.iterdir():
    if p.suffix in ['.json','.txt','.md','.tsv']:
        t=p.read_text(encoding='utf-8')
        if p.suffix=='.txt':t='\n'.join(line.rstrip(' ') for line in t.splitlines())+'\n'
        p.write_text(t,encoding='utf-8',newline='\n')
print(json.dumps(v,ensure_ascii=False))
