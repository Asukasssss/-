"""Presentation only: immutable aggregate input, no new patient statistics."""
from pathlib import Path
import io,json,hashlib,subprocess,csv
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import fitz
ROOT=Path(__file__).resolve().parents[1]
SHA='8952c54b61709c875ec83881b66ee43bd1de59f3'
BASE='results/BRCA/07_INTEGRATION/20260922T140000Z_report_v1/appendix/'
OUT=ROOT/'results/BRCA/07_INTEGRATION/20260922T151000Z_AASS_presentation_v1'
assert not OUT.exists(), 'Refuse overwrite'
OUT.mkdir(parents=True)
manifest=[]
def raw(path):
 b=subprocess.check_output(['git','show',SHA+':'+path],cwd=ROOT);manifest.append(dict(path=path,commit=SHA,sha256=hashlib.sha256(b).hexdigest()));return b
def table(name):
 d=pd.read_csv(io.BytesIO(raw(BASE+name+'.tsv')),sep='\t',dtype=str,keep_default_na=False)
 d=d[d.metabolite_key.eq('KEGG:C00047')] if name=='metabolites' else d[d.gene.eq('AASS')]
 if 'partition' in d:d=d[d.partition.eq('ALL')]
 d.to_csv(OUT/(name+'.tsv'),sep='\t',index=False);return d
assoc=table('association');rna=table('rna');met=table('metabolites');ss=table('source_status');sub=table('subtypes')
profiles=[table('Wu2021_old'),table('Pal2021_reprocessed_old')]
font_manager.fontManager.addfont('C:/Windows/Fonts/msyh.ttc');plt.rcParams.update({'font.family':'Microsoft YaHei','axes.unicode_minus':False,'font.size':11,'pdf.fonttype':42})
fig,axes=plt.subplots(2,2,figsize=(15,11));fig.subplots_adjust(left=.13,right=.97,top=.87,bottom=.15,wspace=.55,hspace=.55)
fig.suptitle('AASS—赖氨酸｜BRCA现有结果',fontsize=23,y=.97,fontweight='bold')
fig.text(.5,.925,'配对变化、肿瘤内部关联与细胞来源是三个不同问题；本图不新增统计检验',ha='center',color='#52616B')
ax=axes[0,0];ax.set_title('A  同患者肿瘤—正常：多数降低',loc='left',fontweight='bold')
for y,(label,low,high) in enumerate([('赖氨酸',34,11),('AASS RNA',37,8)]):
 ax.barh(y,low,color='#087F8C');ax.barh(y,high,left=low,color='#C8D6DF');ax.text(low/2,y,f'{low}/45 降低',ha='center',va='center',color='white');ax.text(low+high/2,y,f'{high}升高',ha='center',va='center',fontsize=9)
ax.set_yticks([0,1],['赖氨酸','AASS RNA']);ax.set_xlim(0,45);ax.set_xlabel('作者明确配对病例数；两行不是同一测量尺度');ax.invert_yaxis()
ax.text(0,-.35,'赖氨酸 q=0.000539；RNA q=2.10e-8',transform=ax.transAxes,fontsize=11)
ax=axes[0,1];ax.set_title('B  肿瘤内部：正相关保留',loc='left',fontweight='bold')
for i,(_,r) in enumerate(assoc.iterrows()):
 e,l,u=map(float,[r.effect,r.ci_lower,r.ci_upper]);ax.errorbar(e,i,xerr=[[e-l],[u-e]],fmt='o',capsize=4,color='#087F8C');ax.text(.68,i,f'n=60\nP={float(r.p_value):.4f}\nq={float(r.q_value):.4f}',va='center',fontsize=10)
ax.set_yticks([0,1],['主分析','可用值分析']);ax.set_xlim(-.1,.95);ax.set_ylim(1.6,-.6);ax.axvline(0,color='#ADB5BD',ls='--');ax.set_xlabel('Spearman ρ及原95%点区间（不是倍数）')
ax.text(0,-.35,'两次ρ均=0.395；q跨0.05不代表效应消失',transform=ax.transAxes,fontsize=11)
labels={'Fibroblasts':'成纤维','Nonmalignant_epithelial':'非恶性上皮','Malignant_epithelial':'恶性上皮','Perivascular':'血管周','Endothelial':'内皮','B_cells':'B细胞','T_cells':'T细胞','T_NK_cells':'T/NK细胞','T_NK':'T/NK','Myeloid':'髓系','Mast_cells':'肥大细胞','Plasma_cells':'浆细胞','Unresolved_stromal':'未定基质'}
for j,d in enumerate(profiles):
 ax=axes[1,j];d=d[d.status.eq('DONE')].copy();d['num']=d.effect.astype(float);d=d.sort_values('num');cohort=['Wu','Pal再注释'][j]
 ax.barh(range(len(d)),d.num,color=['#087F8C' if x=='Fibroblasts' else '#BDCCD4' for x in d.celltype]);ax.set_yticks(range(len(d)),[labels.get(x,x) for x in d.celltype],fontsize=10);ax.set_xlabel('供者等权平均 log1p(counts/库大小×10⁴)',fontsize=10)
 r=ss.iloc[j];ax.set_title(f'{"CD"[j]}  {cohort}：成纤维平均最高',loc='left',fontweight='bold');ax.text(0,-.30,f'成纤维：{r.top_n}个来源供者标签；检出比例 {float(r.top_detection):.1%}\n原最高类别重抽样保持率 {float(r.bootstrap):.1%}',transform=ax.transAxes,fontsize=10)
for ax in axes.flat:ax.spines[['top','right']].set_visible(False)
for ext in ['png','pdf','svg']:fig.savefig(OUT/('AASS_summary.'+ext),dpi=180)
plt.close(fig)
path='results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/all156_expression_celltype_pairs.pdf'
src=fitz.open(stream=raw(path),filetype='pdf');index=pd.read_csv(io.BytesIO(raw(path.rsplit('/',1)[0]+'/gene_page_index.tsv')),sep='\t');page=int(index.loc[index.gene.eq('AASS'),'page'].iloc[0]);assert page==13
umap=fitz.open();umap.insert_pdf(src,from_page=page-1,to_page=page-1);umap.save(OUT/'AASS_UMAP.pdf');umap[0].get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(OUT/'AASS_UMAP.png')
report=fitz.open(OUT/'AASS_summary.pdf');report.insert_pdf(umap);report.save(OUT/'AASS_图册.pdf');report.close();umap.close();src.close()
text='''# AASS—赖氨酸：当前BRCA结果

## 问题与输入
是否同时存在代谢物变化、AASS表达变化、患者间联系及可解释细胞来源？固定输入8952c54（其科学统计冻结于8663004）；无新增P/q，无逐患者/细胞数据导出。

## 实际结果
|层次|结果|含义|
|---|---|---|
|赖氨酸配对|34/45降低（75.6%）；平均处理值差−0.383；P=0.0001424，q=0.0005392|肿瘤相对自己正常组织通常较低；不是浓度降低38.3%|
|赖氨酸可用值|仍45对；P同上，q=0.0009767|该项无配对因可用值限制被排除；不同BH族使q不同|
|AASS配对RNA|37/45降低（82.2%）；平均作者尺度差−0.929，95%区间−1.169至−0.688；q=2.10e-8|RNA降低，不能推断蛋白、酶活或抑癌作用|
|肿瘤内部关联|60个作者病例；rho=0.394554，95%区间0.1445至0.5984；P=0.0025，q=0.04554|肿瘤之间AASS较高者，赖氨酸也倾向较高|
|关联可用值|同60个；rho相同；P=0.0024，q=0.051|效应未变；置换随机种子与全族P分布不同，不把阈值跨越解释为失效|
|Wu来源|成纤维最高；23来源标签；平均检出8.61%；排名保持88.8%|描述表达背景，不是独有作用细胞|
|Pal再注释来源|成纤维最高；30来源标签；平均检出5.06%；排名保持72.6%|未达到两研究均≥80%的稳定规则；整体检出不高|

## 图怎么读
AASS_summary：配对人数、关联区间、两研究全可评估细胞类平均表达。C/D各自横轴，不能跨研究比较绝对柱长。UMAP直接提取已有图册第13页；作者/重注释来源与颜色尺度沿用原图，不重算坐标，不作富集显著性证明。

## 分型背景
Wu ER+与HER2+自身可评估类别中成纤维最高；TNBC非恶性上皮略高于成纤维，差约0.00159。限制三型共有的7个类别后三型均成纤维最高。两种比较范围不同，不称已经证明亚型差异。

## 我的推理与限制
值得保留的是“赖氨酸和AASS RNA在疾病比较中均较低，肿瘤内部又存在正关联，并有成纤维表达背景”这一组合。疾病组间比较与肿瘤内部相关不矛盾。
不能据此说AASS使赖氨酸升高/降低、成纤维细胞主导了代谢变化，或应该激活/抑制AASS。细胞数及状态、共同调控等解释尚未区分；没有本轮专属组成调整或外部复现。
目前定位：有明确患者线索、成纤维背景值得讨论的候选，尚非机制闭合或治疗靶点。保留全候选池，不据本图给AASS唯一第一名。

## 复现与当前决定
`python code/brca_aass_presentation_v1.py`；写新固定运行目录，非空/存在即拒绝覆盖。停在现有结果展示，不新增机制分析。
'''
(OUT/'README_CN.md').write_text(text,encoding='utf-8')
(OUT/'figure_caption_CN.md').write_text('AASS现有汇总的描述性展示。A为方向人数，B为已有Spearman效应与点区间，C/D为各研究供者等权平均表达，横轴分开。全量源行见各TSV。UMAP为固定图册第13页直接复用。所有统计版本/检验族见源表，不新增P/q。',encoding='utf-8')
pd.DataFrame(manifest).to_csv(OUT/'source_manifest.tsv',sep='\t',index=False)
(OUT/'analysis_spec.json').write_text(json.dumps(dict(input_commit=SHA,gene='AASS',relation='KEGG:C00047|AASS',new_tests=0,new_embedding=False,source_page=13),indent=2),encoding='utf-8')
(OUT/'validation.json').write_text(json.dumps(dict(source_relation_rows=len(assoc),paired_rna_rows=len(rna),source_studies=len(ss),new_P_q=0,visual_review='PENDING'),indent=2),encoding='utf-8')
pd.DataFrame([dict(path=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(OUT.iterdir())]).to_csv(OUT/'checksums.tsv',sep='\t',index=False)
print(OUT)
