"""Report every predeclared parameter setting; never select a winning model."""
from pathlib import Path
import sys,json,hashlib
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.backends.backend_pdf import PdfPages
R=Path(sys.argv[1]);F=R/'figures';F.mkdir(exist_ok=True);font=Path('C:/Windows/Fonts/msyh.ttc')
if font.exists():plt.rcParams['font.family']=FontProperties(fname=str(font)).get_name()
plt.rcParams.update({'axes.unicode_minus':False,'pdf.fonttype':42,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
grid=pd.read_csv(R/'parameter_grid.tsv',sep='\t');configs=grid.config.tolist()
q=pd.concat([pd.read_csv(R/(c+'_topic20.tsv'),sep='\t') for c in configs],ignore_index=True)
s=pd.concat([pd.read_csv(R/(c+'_summary.tsv'),sep='\t') for c in configs],ignore_index=True)
allgenes=pd.concat([pd.read_csv(R/(c+'_all_genes.tsv'),sep='\t') for c in configs],ignore_index=True)
assert len(q)==160 and not q.duplicated(['config','gene']).any()
assert len(allgenes)==48000 and not allgenes.duplicated(['config','gene']).any()
assert bool(s.loc[s.config.eq('baseline'),'baseline_exact_labels'].iloc[0])
ok=q.status.eq('DONE');q['top10_percent']=ok & q.rank_fraction.le(.1);q['top20_percent']=ok & q.rank_fraction.le(.2)
summary=[]
for gene,a in q.groupby('gene'):
    b=a[a.config=='baseline'].iloc[0];valid=a[a.status=='DONE'];n=len(valid)
    summary.append(dict(gene=gene,configs_planned=8,configs_evaluable=n,baseline_module=b.module,baseline_size=b.module_size,baseline_rank=b['rank'],baseline_rank_fraction=b.rank_fraction,best_observed_rank_fraction=valid.rank_fraction.min(),worst_observed_rank_fraction=valid.rank_fraction.max(),median_rank_fraction=valid.rank_fraction.median(),top10_configs=int(a.top10_percent.sum()),top20_configs=int(a.top20_percent.sum()),min_kME=valid.kME.min(),max_kME=valid.kME.max(),min_Pal_kME=valid.Pal_kME.min(),max_Pal_kME=valid.Pal_kME.max(),Pal_fixed_module_top10_configs=int(valid.Pal_rank_fraction.le(.1).sum()),minimum_module_jaccard=valid.baseline_module_jaccard.min(),status='DONE' if n else 'NOT_EVALUABLE'))
summ=pd.DataFrame(summary);summ.to_csv(R/'topic20_parameter_summary.tsv',sep='\t',index=False);q.to_csv(R/'topic20_all_parameters.tsv',sep='\t',index=False);s.to_csv(R/'network_parameter_summary.tsv',sep='\t',index=False)
gene_order=['LYPLA1','LYPLA2','PCYT2','ETNK1','GPCPD1']+[g for g in sorted(q.gene.unique()) if g not in ['LYPLA1','LYPLA2','PCYT2','ETNK1','GPCPD1']]
mat=q.pivot(index='gene',columns='config',values='rank_fraction').reindex(index=gene_order,columns=configs)
fig,ax=plt.subplots(figsize=(13,9));cmap=plt.get_cmap('viridis_r').copy();cmap.set_bad('#eeeeee');im=ax.imshow(mat.to_numpy(),vmin=0,vmax=1,cmap=cmap,aspect='auto')
ax.set_xticks(range(8));ax.set_xticklabels(configs,rotation=35,ha='right');ax.set_yticks(range(20));ax.set_yticklabels(gene_order)
for i in range(20):
    for j in range(8):
        v=mat.iloc[i,j];rgb=cmap(0 if pd.isna(v) else v)[:3];lum=sum(a*b for a,b in zip(rgb,[.2126,.7152,.0722]));ax.text(j,i,'×' if pd.isna(v) else f'{100*v:.0f}%',ha='center',va='center',fontsize=9,color='#666666' if pd.isna(v) else ('white' if lum<.5 else 'black'))
fig.colorbar(im,ax=ax,label='模块内排名 / 模块大小（越小越靠前）');ax.set_title('专题20基因：全部8组预设参数的模块内排名',pad=15)
fig.text(.03,.015,'×：未进入冻结输入或未归入有效模块；8组参数不是8次独立验证。',fontsize=10);fig.tight_layout(rect=(0,.04,1,1));fig.savefig(F/'01_all20_rank_sensitivity.png',dpi=170);fig.savefig(F/'01_all20_rank_sensitivity.pdf')
fig2,axes=plt.subplots(1,3,figsize=(15,5));ly=q[q.gene.eq('LYPLA1')].set_index('config').reindex(configs)
xx=np.arange(8)
axes[0].plot(xx,ly.kME,'o-',label='Wu');axes[0].plot(xx,ly.Pal_kME,'s--',label='Pal 固定模块');axes[0].set(ylabel='kME',title='LYPLA1模块成员相关');axes[0].legend()
axes[1].plot(xx,ly.rank_fraction*100,'o-',color='#277d8e');axes[1].axhline(10,color='grey',ls=':');axes[1].set(ylabel='模块内排名百分比（越低越前）',title='同时记录名次及模块大小',ylim=(0,75))
for i,row in enumerate(ly.itertuples()):
    if pd.notna(row.rank):axes[1].annotate(f'{int(row.rank)}/{int(row.module_size)}',(i,row.rank_fraction*100),xytext=(0,8),textcoords='offset points',ha='center',fontsize=8)
axes[2].plot(xx,ly.baseline_module_jaccard,'o-',color='#a76235');axes[2].set(ylabel='与基线模块的基因集Jaccard',ylim=(-.03,1.08),title='模块内容是否改变')
for a in axes:a.set_xticks(xx);a.set_xticklabels(configs,rotation=60,ha='right',fontsize=8)
fig2.suptitle('LYPLA1：参数敏感性不是以排名选最优模型');fig2.tight_layout(rect=(0,0,1,.93));fig2.savefig(F/'02_LYPLA1_parameter_details.png',dpi=170);fig2.savefig(F/'02_LYPLA1_parameter_details.pdf')
with PdfPages(R/'BRCA_WGCNA_parameter_sensitivity.pdf') as pdf:pdf.savefig(fig);pdf.savefig(fig2)
plt.close('all')
with pd.ExcelWriter(R/'BRCA_WGCNA_parameter_sensitivity.xlsx',engine='openpyxl') as w:
    for name,d in [('参数',grid),('20基因概括',summ),('20基因全部方案',q),('LYPLA1各方案',ly.reset_index()),('网络概况',s),('全基因全部方案',allgenes)]:d.to_excel(w,sheet_name=name,index=False)
    for ws in w.book.worksheets:
        ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions
        for col in ws.columns:ws.column_dimensions[col[0].column_letter].width=min(40,max(16,len(str(col[0].value))+1))
ls=summ[summ.gene.eq('LYPLA1')].iloc[0]
rows='\n'.join(f"|{r.Index}|{int(r.module_size)}|{int(r.rank)}|{r.rank_fraction:.1%}|{r.kME:.3f}|{r.Pal_kME:.3f}|{r.Pal_rank_fraction:.1%}|{r.baseline_module_jaccard:.3f}|" for r in ly.itertuples())
retained=summ[summ.configs_evaluable.gt(0)].sort_values(['top10_configs','top20_configs','median_rank_fraction'],ascending=[False,False,True])
other='\n'.join(f"|{r.gene}|{r.configs_evaluable}|{r.top10_configs}|{r.top20_configs}|{r.best_observed_rank_fraction:.1%}～{r.worst_observed_rank_fraction:.1%}|" for r in retained.itertuples())
text=f'''# BRCA WGCNA参数敏感性 v1

## 本轮问题

用户要求尝试不同WGCNA参数，观察感兴趣基因是否进入模块前列。本轮做有范围的敏感性比较，保留所有设置，不按LYPLA1的最好排名决定主模型。开始前固定8个设置，见parameter_grid.tsv与analysis_spec.json。

## 输入与范围

复用20260923T043651Z_malignant_wgcna_v1的Wu 20位供者、160个不重叠metacells、6000基因；Pal固定重叠5992基因、32位供者256个metacells。没有重新抽细胞、放宽基因入网、添加正常细胞或改变归一化。全部专题20基因都记录；不在固定输入的对象保留NOT_EVALUABLE。

在基线power12/deepSplit2/merge0.25/minSize30基础上，每次只改变一个参数：power10/14、deepSplit1/3、merge0.15/0.35、minSize60。共8组，不是全因子网格。相同power使用同一TOM缓存，减少重复计算；缓存和网络对象只在服务器私有目录。

## 实际结果

LYPLA1在{ls.configs_evaluable}组可评估方案中，进入模块前10%的方案数为{ls.top10_configs}，前20%为{ls.top20_configs}。观察到的排名比例范围为{ls.best_observed_rank_fraction:.1%}～{ls.worst_observed_rank_fraction:.1%}。最好名次是多方案中选出的描述性值，不是独立验证或新的显著性。

|设置|模块大小|LYPLA1名次|名次/大小|Wu kME|Pal kME|Pal固定模块排名比例|相对基线Jaccard|
|---|---|---|---|---|---|---|---|
{rows}

全部有可评估结果的专题基因：

|基因|可评估方案|前10%方案数|前20%方案数|排名比例范围|
|---|---|---|---|---|
{other}

## 新手解释

模块内第10名在50基因模块中属于前20%，在500基因模块中属于前2%，所以不能只比较名次。本轮同时列出kME、模块大小和名次比例。kME为基因与模块概括表达的相关，不是表达倍数或治疗作用。

Jaccard为该基因所在模块与基线模块的基因集合重叠/并集。若名次改善但Jaccard很低，说明换成了不同模块内容，不能说同一个模块内稳定升到了前列。不同设置的颜色不是固定生物学身份。

Pal kME使用该方案在Wu得到的模块成员集合，不是重新在Pal发现同一模块；也不是各方案的正式模块保留检验。正式100次置换保留仅在基线批次计算。

LYPLA1在Pal固定模块中的前10%方案数为{ls.Pal_fixed_module_top10_configs}/8。该相对排名受其余模块成员在Pal中的表现影响；不能只因排名较高就忽略约0.5的kME或宣称两研究都把它独立识别为hub。

## 限制与反证

- 8组参数使用同一输入，不是8个独立实验。前10%/20%为描述性标记，不赋予靶点资格。
- 范围不涵盖全部参数、metacell抽样、基因集与归一化，因此不能证明对所有设置都稳定。
- 不因目标排名挑选power；10/14是围绕基线的敏感性设置，不宣称都满足相同拓扑筛选条件。
- 未对每个新模块重新做GO或全网供者留出；不能把基线模块的功能富集直接套到改变后的模块。
- 缺测/过滤不是生物学阴性。grey不是有效功能模块，完整记录保留。
- 全部6000基因×8组成员表均交付；本轮不新增P/q、不改旧CAMP统计。

## 当前决定

主模型继续保留原基线。根据完整参数范围评价“是否经常居前”，而不是展示最有利的一次。即使有方案进入前列，也仍需结合外部表现及原患者证据；不能称为WGCNA证明的驱动或可用治疗靶点。

## 下一步

本轮到完整参数比较交付为止，不继续扩大搜索直到某个基因成为hub。需要改变细胞范围或研究问题时另立分析版本。

## 复现与文件

在server165新建独占RUN，并复制parameter_grid.tsv到RUN/public。使用基线Rlib与private/preservation_input.rds：

```text
Rscript brca_wgcna_parameter_grid_v1.R RUN BASELINE_RUN 10
Rscript brca_wgcna_parameter_grid_v1.R RUN BASELINE_RUN 12
Rscript brca_wgcna_parameter_grid_v1.R RUN BASELINE_RUN 14
```

仅下载public/汇总文件，运行`python brca_wgcna_grid_present_v1.py RESULT_DIRECTORY`。PDF、Excel与两个PNG为展示入口。统计列不含P/q：effect=kME或rank_fraction，unit=metacell共表达描述；供者/细胞源表不公开。
'''
(R/'README_CN.md').write_text(text,encoding='utf-8')
(R/'validation.json').write_text(json.dumps(dict(status='PASS',all_8_configs=True,topic_rows=160,all_gene_rows=48000,unique_config_gene=True,baseline_reproduced_exact_module_labels=True,no_new_P_q=True,ranking_not_model_selection=True),indent=2))
print(summ[summ.gene.isin(['LYPLA1','LYPLA2','PCYT2','ETNK1','GPCPD1'])].to_string(index=False));print(ly.to_string())
