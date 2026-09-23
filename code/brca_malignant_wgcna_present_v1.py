"""Aggregate-only presentation/validation; no source matrices required."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties

R=Path(sys.argv[1]); F=R/'figures';F.mkdir(exist_ok=True)
font=Path('C:/Windows/Fonts/msyh.ttc')
if font.exists():plt.rcParams['font.family']=FontProperties(fname=str(font)).get_name()
plt.rcParams.update({'axes.unicode_minus':False,'pdf.fonttype':42,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
def read(n):return pd.read_csv(R/n,sep='\t')
m=read('gene_modules_external_kME.tsv');mods=read('module_summary.tsv');cov=read('coverage.tsv');soft=read('soft_power.tsv')
ly=m[m.gene.eq('LYPLA1')]; assert len(ly)==1
target=ly.iloc[0]['module']; members=m[m.module.eq(target)].sort_values('kME',ascending=False)
en=read('GO_BP_all_modules.tsv'); own=en[en.module.eq(target)].sort_values(['q_all_module_terms','p'])
pres_done=(R/'preservation_Z.tsv').exists() and (R/'preservation_spec.json').exists()
z=read('preservation_Z.tsv') if pres_done else pd.DataFrame(columns=['module','Zsummary.pres'])
loo=read('LYPLA1_leave_one_donor_out.tsv') if (R/'LYPLA1_leave_one_donor_out.tsv').exists() else pd.DataFrame()
summary=json.loads((R/'network_summary.json').read_text());ps=json.loads((R/'preservation_spec.json').read_text()) if pres_done else dict(status='PARTIAL',reason='100-permutation preservation output pending;not used for conclusions',reference_donors=int(cov.iloc[0].eligible_donors),external_donors=int(cov.iloc[1].eligible_donors))
topic_genes='ABHD12 CEPT1 CHKA CHKB CHPT1 ENPP2 ETNK1 ETNK2 GDPD5 GPCPD1 LPCAT1 LPCAT2 LPCAT3 LPCAT4 LYPLA1 LYPLA2 PCYT1A PCYT1B PCYT2 SLC22A2'.split()
topic=pd.DataFrame({'gene':topic_genes}).merge(m,on='gene',how='left').merge(read('Wu_gene_filter.tsv')[['gene','eligible','detection','donors_detected','network_selected']],on='gene',how='left')
topic['status']=np.where(topic.module.notna(),'DONE','NOT_EVALUABLE');topic['reason']=np.where(topic.module.notna(),'','did not enter network under common gene filtering/ranking;not biological negative')
topic.to_csv(R/'choline_ethanolamine_20_module_lookup.tsv',sep='\t',index=False)
pages=[]
def save(fig,name):
    fig.savefig(F/(name+'.png'),dpi=170,bbox_inches='tight',facecolor='white')
    fig.savefig(F/(name+'.pdf'),bbox_inches='tight',facecolor='white');pages.append(fig)

fig,ax=plt.subplots(1,2,figsize=(12,5));fig.suptitle('恶性上皮共表达网络：全部模块')
ax[0].plot(soft.Power,soft['SFT.R.sq'],'o-',color='#176b88');ax[0].axhline(.8,ls='--',color='grey');ax[0].axvline(summary['power'],ls=':',color='#bb4b27');ax[0].set(xlabel='软阈值 power',ylabel='Scale-free fit R²',title=f"使用 power={summary['power']}；fallback={summary['power_fallback']}")
mm=mods.sort_values('n_genes');ax[1].barh(mm.module,mm.n_genes,color=['#d07138' if a==target else '#7ea5b7' for a in mm.module]);ax[1].set(xlabel='基因数',title=f'LYPLA1所在模块：{target}')
fig.text(.02,.01,'模块颜色是标签；不代表功能。未归入模块的grey单独保留。',fontsize=9);fig.tight_layout(rect=(0,.04,1,.94));save(fig,'01_network_overview')

if target!='grey':
    top=members.head(20).copy()
    if 'LYPLA1' not in set(top.gene):top=pd.concat([top,members[members.gene.eq('LYPLA1')]])
    top=top.sort_values('kME');fig,ax=plt.subplots(1,2,figsize=(12,7));fig.suptitle(f'LYPLA1模块：{target}，{len(members)}个基因')
    ax[0].barh(top.gene,top.kME,color=['#d07138' if a=='LYPLA1' else '#43849e' for a in top.gene]);ax[0].set(xlabel='Wu kME（基因与模块概括表达的相关）',title='模块内kME最高20个，另保留LYPLA1')
    ax[1].scatter(members.kME,members.Pal_kME,s=14,alpha=.5,color='#43849e');ax[1].scatter(ly.kME,ly.Pal_kME,color='#d07138',s=75);ax[1].annotate('LYPLA1',(ly.kME.iloc[0],ly.Pal_kME.iloc[0]),xytext=(5,7),textcoords='offset points');ax[1].set(xlabel='Wu kME',ylabel='Pal 固定模块基因集kME',title='跨研究描述性比较')
    fig.text(.02,.01,'kME高表示共表达中心性较高，不等于上游调控、蛋白互作或治疗依赖。',fontsize=9);fig.tight_layout(rect=(0,.04,1,.94));save(fig,'02_LYPLA1_module')
    import textwrap
    show=own.head(15).iloc[::-1].copy();fig,ax=plt.subplots(figsize=(12,7));fig.suptitle('LYPLA1模块的GO生物过程：按全模块联合q排序')
    labels=['\n'.join(textwrap.wrap(t,55)) for t in show.description];q=show.q_all_module_terms.clip(lower=1e-300)
    ax.scatter(-np.log10(q),range(len(show)),s=20+show.overlap*3,c=np.where(q<.05,'#43849e','#aaaaaa'));ax.set_yticks(range(len(show)));ax.set_yticklabels(labels,fontsize=9);ax.axvline(-np.log10(.05),color='grey',ls='--');ax.set(xlabel='−log10（所有模块×GO条目联合BH q）',title='点面积随重叠基因数增加；完整结果含未显著条目')
    fig.tight_layout(rect=(0,.02,1,.94));save(fig,'03_LYPLA1_GO')

fig,ax=plt.subplots(1,2,figsize=(12,5));fig.suptitle('外部保留与供者敏感性')
zc='Zsummary.pres';sz=z[~z.module.isin(['grey','gold'])].sort_values(zc)
ax[0].barh(sz.module,sz[zc],color=['#d07138' if a==target else '#7ea5b7' for a in sz.module]);ax[0].axvline(2,color='grey',ls=':');ax[0].axvline(10,color='grey',ls='--');ax[0].set(xlabel='Zsummary（探索性保留指标）',title=f"供者汇总：Wu {ps['reference_donors']} / Pal {ps['external_donors']}；100次置换")
if not pres_done:ax[0].text(.5,.5,'置换保留结果尚未完成\n本轮不据此下结论',ha='center',va='center',transform=ax[0].transAxes)
if len(loo):
    ax[1].plot(range(1,len(loo)+1),loo.kME_excluding_self,'o',color='#43849e');ax[1].set(xlabel='留出一次（不展示供者身份）',ylabel='LYPLA1与模块其余基因概括表达的kME',title='每次留出一位Wu供者；固定模块基因集')
else:ax[1].text(.1,.5,'LYPLA1未进入有效模块：不可评估')
fig.text(.02,.01,'留出检查未重新聚类，不能解释为模块归属稳定性；Z受模块大小影响，不是靶点验证。',fontsize=9);fig.tight_layout(rect=(0,.05,1,.93));save(fig,'04_preservation_sensitivity')
with PdfPages(R/'BRCA_LYPLA1_WGCNA_figures.pdf') as pdf:
    for fig in pages:pdf.savefig(fig,bbox_inches='tight')
for fig in pages:plt.close(fig)
with pd.ExcelWriter(R/'BRCA_LYPLA1_WGCNA_results.xlsx',engine='openpyxl') as w:
    cov.to_excel(w,sheet_name='覆盖',index=False);mods.to_excel(w,sheet_name='全部模块',index=False);m.to_excel(w,sheet_name='全部基因',index=False);members.to_excel(w,sheet_name='LYPLA1模块',index=False);own.to_excel(w,sheet_name='LYPLA1_GO全量',index=False);z.to_excel(w,sheet_name='供者层面外部保留',index=False);loo.to_excel(w,sheet_name='固定模块留出',index=False);soft.to_excel(w,sheet_name='软阈值',index=False)
    topic.to_excel(w,sheet_name='专题20基因归属',index=False)
    for ws in w.book.worksheets:
        ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions
        for col in ws.columns:ws.column_dimensions[col[0].column_letter].width=min(55,max(15,len(str(col[0].value))+2))
# Independently check the R BH vector; not a new testing family.
p=en.p.to_numpy();ix=np.argsort(p,kind='stable');expected=np.empty(len(p));expected[ix]=np.minimum(1,np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);err=float(np.max(np.abs(expected-en.q_all_module_terms)))
assert err<1e-10;assert m.gene.is_unique;assert len(m)==summary['network_genes'];assert int(mods.n_genes.sum())==len(m);assert cov.cell_overlap.eq(0).all();assert (cov.metacells==cov.eligible_donors*8).all()
validation=dict(status='PASS',delivery_status='DONE' if pres_done else 'PARTIAL',preservation_completed=pres_done,gene_unique=True,module_counts_sum=int(mods.n_genes.sum()),BH_max_absolute_error=err,zero_metacell_cell_overlap=True,equal_metacells_per_donor=8,source_matrices_recomputed_locally=False,full_network_leave_one_out='NOT_RUN',causal_validation='NOT_RUN')
(R/'validation.json').write_text(json.dumps(validation,indent=2),encoding='utf-8')
print(json.dumps({'LYPLA1':ly.to_dict('records'),'module_size':len(members),'GO_q05':int((own.q_all_module_terms<.05).sum()),'top_GO':own[['description','q_all_module_terms']].head(8).to_dict('records'),'preservation':z[z.module.eq(target)].to_dict('records'),'validation':validation},ensure_ascii=False,indent=2))
run=R.name
standard=pd.DataFrame(index=en.index)
for key,value in dict(cancer='BRCA',cohort='Wu2021',stage_id='06_EXTERNAL',run_id=run,analysis_version='malignant_wgcna_v1',analysis_type='module_GO_BP_ORA',metabolite_key='NA',metabolite_name='NA',gene='NA',unit='annotated_network_gene',n=en.module_annotated_genes,n_reference=en.universe_genes,effect_type='overlap_gene_count',effect=en.overlap,ci_lower=np.nan,ci_upper=np.nan,p_value=en.p,q_value=en.q_all_module_terms,test_family='all_non_grey_modules_x_GO_BP_terms',family_n_evaluable=len(en),status='DONE',reason='ORA;not causal;not patient hypothesis test',source_id='org.Hs.eg.db_GO.db').items():standard[key]=value
standard['module']=en.module;standard['GO']=en.GO;standard.to_csv(R/'results.tsv',sep='\t',index=False,na_rep='NA')
sig=own[own.q_all_module_terms<.05]
lipid=own[own.description.str.contains('lipid|phospholipid|glycerophosph|acyl|cholesterol|fatty acid',case=False,regex=True)]
lipid.to_csv(R/'LYPLA1_module_lipid_keyword_index.tsv',sep='\t',index=False)
spec=json.loads((R/'analysis_spec.json').read_text());spec['external_gene_handling']='Wu network retains all eligible genes; Pal uses uniquely matching symbols only; missing symbols listed without alias guessing';spec['external_preservation']=ps;spec['software']=summary
(R/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf-8')
zly=z[z.module.eq(target)];zval=float(zly[zc].iloc[0]) if len(zly) else np.nan
loostr=f"{loo.kME_excluding_self.min():.3f}～{loo.kME_excluding_self.max():.3f}" if len(loo) else '不可评估'
go_lines='\n'.join(f"|{row.description}|{row.overlap}|{row.q_all_module_terms:.3g}|" for row in own.head(10).itertuples())
lipid_lines='\n'.join(f"- {row.description}：重叠{row.overlap}基因，联合q={row.q_all_module_terms:.3g}。" for row in lipid.head(8).itertuples())
text=f'''# BRCA恶性上皮共表达：LYPLA1模块 v1

交付状态：{'DONE：100次置换保留已完成。' if pres_done else 'PARTIAL：构网、GO、固定模块Pal kME和供者留出已完成；100次置换保留仍在服务器计算，本版不含其结论。'}

## 本轮问题

在作者标注的恶性上皮内构建共表达模块，查明LYPLA1的模块归属、同模块基因、功能富集和跨研究表现。研究目的为潜在靶点的数据挖掘，不是因果或治疗依赖验证。

## 输入与范围

- 发现：Wu2021原作者celltype_major；外部：Pal2021_reprocessed，使用Chen2026对Pal计数的重注释，不冒充Pal原作者标签。均复用冻结标签，未重聚类或重新判定恶性。
- 实际覆盖见coverage.tsv；Wu纳入{int(cov.iloc[0].eligible_donors)}个作者供者标签、{int(cov.iloc[0].metacells)}个metacells；Pal纳入{int(cov.iloc[1].eligible_donors)}个作者供者标签、{int(cov.iloc[1].metacells)}个metacells。
- 每位纳入供者至少120个恶性上皮细胞；先最多抽1000细胞用于候选矩阵/PCA，再在供者内部形成8个不重叠的15细胞近邻组。每位供者对网络贡献相同，未筛选LYPLA1高表达细胞。
- Wu按唯一基因符号、非线粒体、候选细胞检出率≥5%、至少3位供者检出、正方差过滤，选log1p归一化方差最高6000基因。LYPLA1自行通过过滤，没有强制入网；不是只用原候选基因建网。
- 候选细胞集每供者上限相同，但实际数可不同，因而基因筛选/PCA不是严格供者等权。最终构网8组/供者等权。供者汇总使用候选细胞集，最多1000/供者。
- 源矩阵、细胞成员、供者表达均留server165。source_manifest.tsv记录源SHA256。原CAMP统计不变。

## 实际结果

|项目|结果|
|---|---|
|构网基因|{len(m)}|
|有效模块（不含grey）|{summary['modules_excluding_grey']}|
|软阈值|{summary['power']}，fallback={summary['power_fallback']}|
|LYPLA1模块|{target}|
|模块基因数|{len(members)}|
|LYPLA1 Wu kME|{ly.kME.iloc[0]:.4f}|
|模块内kME名次|{ly.rank_in_module.iloc[0]}|
|Pal固定模块LYPLA1 kME|{ly.Pal_kME.iloc[0]:.4f}|
|留出一位供者的kME范围（模块排除LYPLA1自身）|{loostr}|
|供者层面外部Zsummary|{f'{zval:.3f}' if pres_done else 'NA：计算未完成'}|
|本模块GO BP联合q＜0.05条目|{len(sig)}|

完整模块和全部基因均交付，不只展示LYPLA1所在模块。模块名称仅为颜色标签。

GO前10项（未显著仍展示）：

|GO生物过程|重叠基因数|全部模块×条目联合q|
|---|---|---|
{go_lines}

脂质关键词索引是方便查阅的附加索引，不改变检验范围或q；只包含名称匹配，不是穷尽的脂质功能判定：

{lipid_lines}

## 新手解释

模块表示RNA经常一起变化的一组基因。kME表示某个基因与该模块概括表达的相关程度；不是表达倍数，也不是干预效应。中心性高只能提示值得进一步比较，不能因此称为驱动靶点。

GO富集表示该模块相较本次构网背景含有更多某类注释基因；不等于在细胞里测到相应代谢物、酶活或通量。若脂质条目没有通过校正，不能因为LYPLA1已知生化功能就把整个模块命名为脂质模块。

外部Zsummary用于探索性共表达结构保留，常用2和10参考线，但它受模块大小、样本覆盖影响，不能视为P值或治疗验证。Pal kME只是固定模块成员的描述性相关，不是重新发现了相同模块。

## 限制与反证

- 使用自定义不重叠metacells加WGCNA {summary['WGCNA']}，未运行hdWGCNA软件包。邻居在同供者内选，PCA分别按研究计算；不存在跨研究混合metacell。
- bicor遇到零MAD的列按预设pearsonFallback=individual回退至Pearson，软件日志中的相应警告被保留。
- 160/256个细胞组不是160/256名独立患者。主网络用于发现；正式外部保留以供者汇总作为输入（Wu {ps['reference_donors']}、Pal {ps['external_donors']}）。作者身份标签未做基因型独立性确认。
- 留出一位供者只重算固定模块的概括表达与LYPLA1相关，没有重建全网络，所以不能称为模块归属重采样稳定率。
- 没有消除亚型、供者、技术批次、细胞周期或CNV的全部影响。不能把跨细胞共表达直接解释为细胞内调控；亚型均值表为描述性，无新亚型P/q。
- Pal无唯一同名记录的基因单列，不拼接别名、不填0；供者保留进一步排除零方差基因。模块保留100次置换，大模块最多抽1000基因，属于探索性精度。
- 无正常上皮参与建网，本轮不重新检验恶性比正常是否高表达。
- GO版本与全部检验族见enrichment_spec.json；BH跨全部模块×全部合格BP条目（包括零重叠）。关键词索引不是单独校正的脂质检验族。
- 已知LYPLA1乳腺癌研究及原患者结果仍是独立证据栏，本轮不改写已有研究的新颖性，也不把共表达补充当作功能实验。

## 当前决定

保留LYPLA1作为已有患者/表达线索上的共表达专题对象。是否与脂质程序相关，以完整GO结果、同模块基因和外部保留共同解释；不为得到某个模块反复改阈值。其他模块及未支持条目完整保留。后续无需自动展开完整机制研究。

## 下一步

先阅读本模块前列基因、GO全表与四张图，判断这条模块信息是否真的增加候选价值。若主要反映普遍翻译/增殖程序，也如实记录，不强行包装成脂质特异性。后续干预或其他队列是独立新问题，不是本轮完成条件。

## 复现命令与记录

在server165新建独占运行目录，复制脚本、冻结的sc117读取helper和Wu/Pal配置；安装WGCNA到运行目录Rlib。配置文件及helper均在仓库code/。依次运行：

```text
Rscript brca_malignant_wgcna_install_v1.R RUN
python3 brca_malignant_metacells_v1.py RUN
Rscript brca_malignant_wgcna_v1.R RUN
Rscript brca_malignant_wgcna_followup_v1.R RUN
```

仅将public/同步到本结果目录，再在有pandas/matplotlib/openpyxl环境运行`python brca_malignant_wgcna_present_v1.py RESULT_DIRECTORY`。网络无需任何LYPLA1专用调参。

首次准备因Pal缺少8个Wu同名基因被严格检查终止；随后明确缺项处理并按原种子和参数重跑，未因结果改变构网参数。错误日志留服务器。当前版本提供完整缺项表。

文件入口：BRCA_LYPLA1_WGCNA_figures.pdf、BRCA_LYPLA1_WGCNA_results.xlsx；全部数值TSV保留原WGCNA字段。共有统计字段的对应：cancer=BRCA、stage=06_EXTERNAL、run={run}；kME是bicor、GO效应为overlap并配hypergeometric P/BH q、preservation是Z统计量，无共同的“重要性P值”。不适用的区间/P/q不补造。

方法参考：[WGCNA](https://doi.org/10.1186/1471-2105-9-559)、[hdWGCNA的metacell方法说明](https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html)。后者只作方法参考，本轮实现有上述明确差异。
'''
(R/'README_CN.md').write_text(text,encoding='utf-8')
