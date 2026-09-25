"""Verify KEGG aggregate results, join lipid pathways, and render reports."""
import argparse,hashlib,json,sys,textwrap
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'runtime/report_dependencies'))
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();out=a.out
def read(name):return pd.read_csv(out/name,sep='\t')
def save(frame,name):frame.to_csv(out/name,sep='\t',index=False,na_rep='NA',lineterminator='\n')
d=read('gene_statistics.tsv');ora=read('kegg_ORA.tsv');g=read('kegg_GSEA.tsv');cov=read('pathway_coverage.tsv');mapping=read('gene_mapping.tsv');members=read('tested_pathway_members.tsv')
spec=json.loads((out/'analysis_spec.json').read_text());v=json.loads((out/'validation.json').read_text())
assert v['status']=='PASS' and len(d)==11058 and d.gene.is_unique
prior=pd.read_csv(ROOT/'results/COAD/06_EXTERNAL/20260925T153907Z_lypla1_fc025_v3/results.tsv',sep='\t').set_index('gene')
for col in ['log2FC','QL_F','p_value','depth_log2FC','depth_p_value']:
    assert np.allclose(d.set_index('gene')[col],prior.loc[d.gene,col],rtol=1e-13,atol=1e-15)
assert np.allclose(hypergeom.sf(ora.overlap-1,ora.background_genes,ora.pathway_tested_genes,ora.selected_genes),ora.p_value,rtol=1e-10,atol=1e-15)
sets=members.groupby('pathway_id').gene.agg(set).to_dict();lookup=d.set_index('gene');universe=set(d.gene)
selected=(d.p_value<.05)&(d.log2FC.abs()>=.25)
assert selected.sum()==434
for _,r in ora.iterrows():
    bg=universe if r.background=='all_tested' else set(mapping.loc[mapping.in_any_pathway,'gene'])
    chosen=set(d.loc[selected&((d.log2FC>0) if r.direction=='higher_in_LYPLA1_high' else (d.log2FC<0)),'gene'])&bg
    assert len(bg)==r.background_genes and len(chosen)==r.selected_genes
    assert (set(str(r.overlap_genes).split(';')) if pd.notna(r.overlap_genes) else set())==sets[r.pathway_id]&chosen
# Independent running-sum reconstruction of every GSEA ES (not its Monte Carlo P).
rank=d[['gene','QL_F','log2FC']].copy();rank['score']=np.sign(rank.log2FC)*np.sqrt(rank.QL_F.clip(lower=0));rank=rank.sort_values(['score','gene'],ascending=[False,True]);weights=np.abs(rank.score.to_numpy())
for _,r in g.iterrows():
    hits=rank.gene.isin(sets[r.pathway_id]).to_numpy();nh=hits.sum();assert nh==r['size']
    step=np.where(hits,weights/max(weights[hits].sum(),1e-300),-1/(len(rank)-nh));running=np.cumsum(step)
    hi=max(0.,running.max());lo=min(0.,running.min());es=hi if hi>-lo else lo
    assert np.isclose(es,r.ES,rtol=1e-8,atol=1e-10),(r.pathway_id,es,r.ES)
    assert set(str(r.leadingEdge).split(';'))<=sets[r.pathway_id]
main=ora[ora.background=='all_tested'];sensitivity=ora[ora.background!='all_tested']
eligible=cov[cov.status=='DONE'];assert len(g)==len(eligible)==v['pathways_evaluable']
links=[]
for _,r in main[main.p_value<.05].iterrows():
    for gene in str(r.overlap_genes).split(';'):
        x=lookup.loc[gene];links.append(dict(method='ORA',pathway_id=r.pathway_id,pathway=r.pathway,direction=r.direction,pathway_P=r.p_value,gene=gene,log2FC=x.log2FC,gene_P=x.p_value,depth_log2FC=x.depth_log2FC,depth_P=x.depth_p_value))
for _,r in g[g.pval<.05].iterrows():
    for gene in str(r.leadingEdge).split(';'):
        x=lookup.loc[gene];links.append(dict(method='GSEA_leading_edge',pathway_id=r.pathway_id,pathway=r.pathway,direction='higher_in_LYPLA1_high' if r.NES>0 else 'lower_in_LYPLA1_high',pathway_P=r.pval,gene=gene,log2FC=x.log2FC,gene_P=x.p_value,depth_log2FC=x.depth_log2FC,depth_P=x.depth_p_value))
save(pd.DataFrame(links),'pathway_gene_links.tsv')
# Official lipid-metabolism category, plus explicitly labelled adjacent pathways.
adjacent={'hsa03320':'PPAR signaling','hsa04975':'Fat digestion and absorption','hsa04979':'Cholesterol metabolism','hsa04976':'Bile secretion','hsa04070':'Phosphatidylinositol signaling','hsa04071':'Sphingolipid signaling','hsa04216':'Ferroptosis'}
focus=cov[cov.lipid_metabolism|cov.pathway_id.isin(adjacent)].copy()
focus['focus_basis']=np.where(focus.lipid_metabolism,'KEGG_LIPID_METABOLISM','RELATED_NOT_LIPID_CATEGORY')
focus=focus.rename(columns={'status':'coverage_status','reason':'coverage_reason'})
for direction,label in [('higher_in_LYPLA1_high','high'),('lower_in_LYPLA1_high','low')]:
    z=main[main.direction==direction][['pathway_id','p_value','overlap_genes','overlap']].rename(columns={c:'ORA_'+label+'_'+c for c in ['p_value','overlap_genes','overlap']})
    focus=focus.merge(z,on='pathway_id',how='left',validate='one_to_one')
    z=sensitivity[sensitivity.direction==direction][['pathway_id','p_value']].rename(columns={'p_value':'ORA_annotated_'+label+'_P'})
    focus=focus.merge(z,on='pathway_id',how='left',validate='one_to_one')
focus=focus.merge(g[['pathway_id','NES','pval','leadingEdge']].rename(columns={'pval':'GSEA_P'}),on='pathway_id',how='left',validate='one_to_one')
focus=focus.sort_values(['focus_basis','pathway_id']);save(focus,'lipid_pathways.tsv')
# Shared statistical schema; q intentionally NA. ORA counts genes, GSEA ranks genes;
# these are not new cell-independent tests or extra patients.
prefix=(ROOT/'templates/statistical_result.tsv').read_text().strip().split('\t');rows=[]
for _,r in ora.iterrows():
    rows.append(dict(analysis_type='KEGG_ORA_'+r.background,effect_type='fold_enrichment',effect=r.fold_enrichment,p_value=r.p_value,pathway_id=r.pathway_id,pathway=r.pathway,direction=r.direction,status=r.status,reason=r.reason,test_family=r.background+'_'+r.direction,family_n_evaluable=len(eligible)))
for _,r in g.iterrows():
    rows.append(dict(analysis_type='KEGG_GSEA',effect_type='NES',effect=r.NES,p_value=r.pval,pathway_id=r.pathway_id,pathway=r.pathway,direction='higher_in_LYPLA1_high' if r.NES>0 else 'lower_in_LYPLA1_high',status=r.status,reason=r.reason,test_family='KEGG_GSEA_all_tested_rank',family_n_evaluable=len(g)))
res=pd.DataFrame(rows)
for k,value in dict(cancer='COAD',cohort='Uhlitz_GSE166555',stage_id='06_EXTERNAL',run_id=out.name,analysis_version=spec['analysis_version'],metabolite_key='NA',metabolite_name='NA',gene='NA',unit='gene_set',n=8,n_reference=8,ci_lower='NA',ci_upper='NA',q_value='NA',source_id='KEGG_hsa_official').items():res[k]=value
save(res[prefix+[c for c in res if c not in prefix]],'results.tsv')
plt.rcParams['font.family']='Microsoft YaHei'
fig,axes=plt.subplots(1,2,figsize=(16,6.8));fig.subplots_adjust(left=.26,right=.98,wspace=1.0,bottom=.12,top=.82)
for ax,direction,title in zip(axes,['higher_in_LYPLA1_high','lower_in_LYPLA1_high'],['高组较高基因','高组较低基因']):
    z=main[(main.direction==direction)&(main.p_value<.05)].nsmallest(8,'p_value').iloc[::-1]
    ax.barh(range(len(z)),-np.log10(z.p_value),color='#c86639' if 'higher' in direction else '#397f9b')
    ax.set_yticks(range(len(z)),[textwrap.fill(x,32) for x in z.pathway],fontsize=9);ax.set_xlabel('−log10(原始P)');ax.set_title(title)
fig.suptitle('KEGG差异基因富集｜各方向P最小的8项\nP<0.05且|log2FC|≥0.25；背景为全部11,058个受检基因',fontsize=14)
fig.savefig(out/'kegg_ORA_top.png',dpi=180,bbox_inches='tight');plt.close(fig)
z=focus[focus.coverage_status=='DONE'].sort_values('GSEA_P',na_position='last');fig,ax=plt.subplots(figsize=(12,max(6,len(z)*.35)))
ax.barh(range(len(z)),z.NES,color=['#c86639' if x>0 else '#397f9b' for x in z.NES]);ax.set_yticks(range(len(z)),[textwrap.fill(x,43) for x in z.pathway],fontsize=9);ax.invert_yaxis();ax.axvline(0,c='gray',lw=.7)
for i,(_,r) in enumerate(z.iterrows()):ax.text(r.NES+(.04 if r.NES>=0 else -.04),i,f'P={r.GSEA_P:.3g}'+(' *' if r.GSEA_P<.05 else ''),ha='left' if r.NES>=0 else 'right',va='center',fontsize=8)
ax.set_xlim(min(-2.5,z.NES.min()-.9),max(2.5,z.NES.max()+.9));ax.set_xlabel('NES：正值偏LYPLA1高组，负值偏低组');ax.set_title('KEGG脂代谢及相邻通路｜GSEA完整可评估条目\n* 原始P<0.05；颜色只表示方向，不表示显著或通路活性',fontsize=13)
fig.tight_layout();fig.savefig(out/'kegg_lipid_GSEA.png',dpi=180,bbox_inches='tight');plt.close(fig)
def pstr(x):return 'NA' if pd.isna(x) else f'{x:.4g}'
ftable='\n'.join(f'|{r.pathway}|{r.tested_gene_count}|{pstr(r.ORA_high_p_value)}|{pstr(r.ORA_low_p_value)}|{pstr(r.NES)}|{pstr(r.GSEA_P)}|' for _,r in focus.iterrows())
txt=f'''# COAD：LYPLA1高低表达的KEGG富集补充

## 本轮问题

按用户要求补做KEGG，沿用原始P<0.05、不使用FDR筛选。与既有Reactome结果并列；更换注释体系不是独立验证，不按哪套P更小选择结论。

## 输入与范围

Uhlitz GSE166555作者CNA上皮，LYPLA1检出阳性细胞按供者内中位数分高低。8位供者、2,528细胞；复用供者配对edgeR结果，11,058受检基因，LYPLA1本身排除。P<0.05且|log₂FC|≥0.25得到434基因：高组较高158、较低276。没有重新分组或计算差异表达。

使用[KEGG官方API](https://www.kegg.jp/kegg/rest/keggapi.html)取得人类通路、基因、成员和分类，来源哈希在source_manifest.tsv。服务器直接下载超时后，官方注释经本机内存转传到服务器，未下载患者矩阵。优先精确主符号匹配，唯一别名可用；歧义与多个受检符号指向同一KEGG ID均排除，去向完整保留。

KEGG人类通路共{v['pathways_total']}项，排除全局/概览图及交集小于15或大于500基因的条目后，{v['pathways_evaluable']}项可评估。疾病命名条目保留，但不视为样本患有该病。{v['genes_annotated']}个受检基因有至少一个KEGG通路注释；上调、下调名单中分别{v['selected_high_annotated']}、{v['selected_low_annotated']}个有通路注释。

ORA主背景与Reactome保持一致，为全部11,058个受检基因；另报告仅KEGG已注释基因的背景敏感性，不择小P。GSEA使用全部受检基因的signed sqrt(QL_F)排序，参数见analysis_spec.json；fgsea默认内部生成的padj被舍弃，不用于筛选或展示。极小负QL_F仅在排序时归零，原效应与P不改动。

## 实际结果

主背景ORA共{int((main.p_value<.05).sum())}个方向—通路组合P<0.05，高组较高{int(((main.p_value<.05)&(main.direction=='higher_in_LYPLA1_high')).sum())}项、高组较低{int(((main.p_value<.05)&(main.direction=='lower_in_LYPLA1_high')).sum())}项。注释背景敏感性有{int((sensitivity.p_value<.05).sum())}项P<0.05。GSEA有{v['GSEA_P05']}项P<0.05，高组端{v['GSEA_high_P05']}项、低组端{v['GSEA_low_P05']}项；两种方法不能相加为独立证据。

脂质相关的主要观察：脂肪消化与吸收（GSEA NES=-1.653，P=0.01229）、PPAR信号（NES=-1.538，P=0.01610）、醚脂代谢（NES=-1.467，P=0.04753）均偏LYPLA1低组。三项对应的差异基因ORA均未达到P<0.05，不能写成两种方法一致验证。脂代谢官方分类中，10项达到15基因门槛，4项覆盖不足；可评估项中仅醚脂代谢GSEA达到本轮原始P阈值。

磷脂酰肌醇信号的高组较低基因ORA P=0.03917，命中ITPKA、ITPKB、ITPR2、PLCD1、PRKCG；改用已注释基因背景后P=0.03788。但该通路GSEA P=0.6098，不能从ORA进一步断言整个通路整体下降。

甘油磷脂代谢（GSEA P=0.6315）、花生四烯酸代谢（P=0.5020）、脂肪酸降解（P=0.3946）未达到P<0.05。脂肪酸生物合成仅13个受检成员、初级胆汁酸生物合成9个、亚油酸代谢10个、α亚麻酸代谢13个，均未达到预先固定的15成员门槛；不为获得结果而临时放宽。

![ORA结果](kegg_ORA_top.png)

脂代谢表按KEGG官方Lipid metabolism分类完整列出，另加明确标记的PPAR、脂肪消化吸收、胆固醇代谢、胆汁分泌、磷脂酰肌醇/鞘脂信号与铁死亡。下表不按P删行；NA表示覆盖规则不满足，不能解释为没有作用。上/下分别表示高组较高/较低基因的ORA；NES正值偏高组、负值偏低组。

|通路|受检基因交集|ORA上P|ORA下P|GSEA NES|GSEA P|
|---|---:|---:|---:|---:|---:|
{ftable}

![脂质相关GSEA](kegg_lipid_GSEA.png)

## 新手解释

ORA检验入选差异基因是否较多落在通路内；GSEA检查通路基因在全排序中的偏向。ORA命中基因与GSEA leading edge含义不同，后者不要求单基因P<0.05。pathway_gene_links.tsv提供所有P<0.05通路的具体基因及其主模型、深度模型数值；lipid_pathways.tsv保留全部脂质条目，包括未过P门槛和不可评估项。

## 限制/反证

全部为原始P探索；共享基因和父子通路不能计作重复支持。分组与总UMI有关：此前高组细胞总UMI中位数9,502，低组13,027.5。434主候选中仅67项在既有深度模型仍同向P<0.05，因此本次通路结果需结合深度敏感性阅读。不能将富集方向直接称作脂质含量、代谢通量或LYPLA1因果调控。CNA为作者基于RNA推断，不是逐细胞DNA证明。注释背景敏感性不是新的供者验证。

## 当前决定

完成KEGG ORA和GSEA，保留完整通路与映射去向；不改变CAMP、Reactome或历史候选名单。只公开汇总，患者/细胞级数据和注释源文件留server165。

## 下一步

按具体通路阅读命中基因及深度敏感性；本批到富集补充收口，不自动扩展机制实验或单细胞模块。

## 复现命令

服务器新独占运行目录中执行prepare_lypla1_kegg_v1.py --out <目录> --commit <代码锁定提交>；采用已核验的官方转传注释时加--cached-sources。随后Rscript run_lypla1_kegg_gsea_v1.R <目录>。本地report_lypla1_kegg_v1.py --out <公开结果目录>重新核对ORA算术、命中集合及全部GSEA累积和效应ES。未独立重算GSEA随机P；参数、代码和来源哈希均保留。
'''
(out/'README_CN.md').write_text(txt,encoding='utf-8',newline='\n')
audit=dict(status='PASS',DE_effects_and_P_unchanged=True,ORA_hypergeometric_and_hits_checked=True,GSEA_all_ES_reconstructed=True,GSEA_P_independently_rerun=False,all_pathway_coverage_retained=True,FDR_selection=False,plots_rendered=True)
(out/'arithmetic_validation.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8',newline='\n')
files=['prepare_lypla1_kegg_v1.py','run_lypla1_kegg_gsea_v1.R','report_lypla1_kegg_v1.py']
(out/'code_manifest.json').write_text(json.dumps([dict(path='code/coad/'+f,sha256=hashlib.sha256((ROOT/'code/coad'/f).read_bytes()).hexdigest()) for f in files],indent=2)+'\n',encoding='utf-8',newline='\n')
for p in out.iterdir():
    if p.suffix in ['.tsv','.json','.txt','.md']:
        t=p.read_text(encoding='utf-8');t='\n'.join(line.rstrip() for line in t.splitlines())+'\n';p.write_text(t,encoding='utf-8',newline='\n')
print(json.dumps(v,ensure_ascii=False));print(json.dumps(audit))
