"""Add existing RNA paired-patient direction counts to the six-cancer gene slide."""
from pathlib import Path
import json,hashlib,shutil,subprocess,io
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'results/BRCA/07_INTEGRATION/20260928T150000Z_ppt_redraw_v1'
PREV=ROOT/'results/BRCA/07_INTEGRATION/20260929T170000Z_page8_genes_six_v2'
OUT=ROOT/'results/BRCA/07_INTEGRATION/20260929T200000Z_page8_patient_consistency_v3'
VIEW=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_page08_patient_consistency_v3')
OUT.mkdir(parents=True,exist_ok=True);VIEW.mkdir(parents=True,exist_ok=True)
sel=pd.read_csv(PREV/'display_genes.tsv',sep='\t');d=pd.read_csv(PREV/'display_cells.tsv',sep='\t')
manifest=pd.read_csv(OLD/'source_manifest.tsv',sep='\t');sources=[];tabs={}
for c in ['brca','coad','pdac']:
 r=manifest[manifest.id.eq(c+'_gene')].iloc[0].to_dict();p=OLD/'sources'/f'{c}_gene.tsv'
 assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'];sources.append(r);tabs[c.upper()]=pd.read_csv(p,sep='\t')
for c,commit,path in [
 ('PRAD','8fe4d512390c5005da36be4920e20ae2517d214a','results/PRAD/03_PATIENT/20260922T142000Z_discovery_v1/paired_RNA_all_families.tsv'),
 ('ccRCC3','02593f745fc966c975ab1a25e3e7ce74f770efc6','results/ccRCC/03_PATIENT/20260925T155500Z_ccrcc3_therapy_v2/rna_results.tsv'),
 ('ccRCC4','02593f745fc966c975ab1a25e3e7ce74f770efc6','results/ccRCC/03_PATIENT/20260925T155000Z_ccrcc4_histology_v2/rna_results.tsv')]:
 b=subprocess.check_output(['git','show',f'{commit}:{path}'],cwd=ROOT);tab=pd.read_csv(io.BytesIO(b),sep='\t');tabs[c]=tab[tab.test_family.eq('RNA_PAIRED_CURRENT')]
 sources.append(dict(id=c+'_rna_primary',git_commit=commit,git_path=path,sha256=hashlib.sha256(b).hexdigest(),bytes=len(b),url=f'https://github.com/Asukasssss/-/blob/{commit}/{path}'))
for name in ['display_cells.tsv','display_genes.tsv']:
 p=PREV/name;sources.append(dict(id='previous_'+name,git_path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size,git_commit='fe76cec'))
fields={
 'BRCA':['v2_paired_RNA_effect','v2_paired_RNA_p_value','v2_paired_RNA_n','v2_paired_RNA_positive_pairs','v2_paired_RNA_negative_pairs','v2_paired_RNA_equal_pairs'],
 'COAD':['RNA_effect','RNA_p_value','RNA_n','RNA_n_up','RNA_n_down','RNA_n_equal'],
 'PDAC':['RNA_effect','RNA_p','RNA_n','RNA_up_pairs','RNA_down_pairs','RNA_equal_pairs']}
records=[]
order=['BRCA','COAD','GBM','PDAC','PRAD','ccRCC3','ccRCC4'];paired=[c for c in order if c!='GBM']
for r in sel.itertuples():
 for c in order:
  a=d[d.gene.eq(r.gene)&d.cohort.eq(c)].iloc[0]
  base=dict(gene=r.gene,cohort=c,main_direction=r.dominant_direction,cohort_effect=a.effect,cohort_p=a.p,same_direction_patients=np.nan,n_pairs=np.nan,up_patients=np.nan,down_patients=np.nan,equal_patients=np.nan,same_fraction=np.nan)
  if a.status!='EVALUABLE':records.append(dict(base,status=a.status));continue
  if c=='GBM':records.append(dict(base,status='NOT_APPLICABLE_UNPAIRED'));continue
  z=tabs[c][tabs[c].gene.eq(r.gene)];assert len(z)==1,(r.gene,c);z=z.iloc[0]
  ec,pc,nc,uc,dc,zc=fields.get(c,['effect','p_value','n','n_up','n_down','n_equal'])
  assert np.isclose(z[ec],a.effect,rtol=0,atol=1e-12) and np.isclose(z[pc],a.p,rtol=1e-12,atol=1e-15) and z[nc]==a.n
  up,down,equal,n=[int(z[k]) for k in [uc,dc,zc,nc]];assert up+down+equal==n
  same=up if r.dominant_direction=='up' else down
  base.update(same_direction_patients=same,n_pairs=n,up_patients=up,down_patients=down,equal_patients=equal,same_fraction=same/n)
  records.append(dict(base,status='DONE'))
counts=pd.DataFrame(records);assert len(counts)==42
counts.to_csv(OUT/'patient_direction_counts.tsv',sep='\t',index=False)
sel.to_csv(OUT/'display_genes.tsv',sep='\t',index=False);d.to_csv(OUT/'display_cells.tsv',sep='\t',index=False)
pd.DataFrame(sources).to_csv(OUT/'source_manifest.tsv',sep='\t',index=False)
shutil.copy2(PREV/'sysu_logo.png',OUT/'sysu_logo.png')

for fp in ['C:/Windows/Fonts/msyh.ttc','C:/Windows/Fonts/msyhbd.ttc']:font_manager.fontManager.addfont(fp)
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
G='#00553B';INK='#183D34';M='#6A7973';UP='#C77462';DOWN='#468B9E';LINE='#DDE7E0'
labels={g:g for g in sel.gene}
f=plt.figure(figsize=(16,9),facecolor='white')
def text(x,y,s,size=15,col=INK,bold=False,ha='left'):return f.text(x,y,s,fontsize=size,color=col,weight='bold' if bold else 'normal',ha=ha,va='center')
def line(x1,x2,y,c=LINE,w=.8):f.add_artist(plt.Line2D([x1,x2],[y,y],transform=f.transFigure,color=c,lw=w))
ax=f.add_axes([.044,.874,.177,.106]);ax.imshow(plt.imread(OUT/'sysu_logo.png'));ax.axis('off')
text(.95,.928,'CAMP  /  组会汇报',10,M,ha='right');line(.05,.95,.876);line(.05,.103,.876,G,2.2)
text(.05,.813,'候选基因的同向变化，在多少患者中重复？',28,G,True)
text(.052,.755,'RNA 肿瘤—参照比较  ·  右侧统计与每行主方向一致的患者数 / 可评估配对数',13,M)
# One common canvas keeps all rows and count columns exactly aligned.
ax=f.add_axes([.05,.236,.90,.472]);ax.set_xlim(0,100);ax.set_ylim(6.15,-1.5);ax.axis('off')
xx=np.linspace(23,56.5,7);cx=63.3;px=np.linspace(72.1,97.6,6)
ax.text(0,-.45,'候选基因',fontsize=11,color=G,weight='bold',va='center')
for x,c in zip(xx,order):ax.text(x,-.45,c,fontsize=10.4,color=G,weight='bold',ha='center',va='center')
ax.text(cx,-.64,'同向显著 / 可评估',fontsize=8.6,color=G,weight='bold',ha='center',va='center');ax.text(cx,-.25,'癌种数',fontsize=9,color=M,ha='center',va='center')
ax.text(84.8,-1.17,'患者一致性',fontsize=13,color=G,weight='bold',ha='center',va='center')
for x,c in zip(px,paired):ax.text(x,-.45,c,fontsize=9.3,color=G,weight='bold',ha='center',va='center')
ax.plot([0,100],[.02,.02],c=LINE,lw=.8)
for x in [59.7,68.4]:ax.plot([x,x],[-1.4,5.8],c=LINE,lw=.8)
for i,r in enumerate(sel.itertuples()):
 y=i+.48;col=UP if r.dominant_direction=='up' else DOWN
 if i%2==0:ax.add_patch(Rectangle((0,y-.44),100,.88,fc='#F5F8F5',ec='none',zorder=-3))
 ax.text(.6,y,labels[r.gene],fontsize=12.3,color=INK,weight='bold',va='center')
 for x,c in zip(xx,order):
  a=d[d.gene.eq(r.gene)&d.cohort.eq(c)].iloc[0]
  if a.status=='EVALUABLE':
   color=UP if a.effect>0 else DOWN;ax.scatter(x,y,s=120,facecolors=color if a.p<.05 else 'white',edgecolors=color,lw=1.4)
  else:ax.text(x,y,'—',fontsize=15,color='#ABB9B1',ha='center',va='center')
 ax.text(cx-2.4,y,'↑' if r.dominant_direction=='up' else '↓',fontsize=16,color=col,ha='center',va='center')
 ax.text(cx+.5,y,f'{r.same_significant_cancers}/{r.evaluable_cancers}',fontsize=16,color=col,ha='center',va='center',weight='bold')
 for x,c in zip(px,paired):
  a=counts[counts.gene.eq(r.gene)&counts.cohort.eq(c)].iloc[0]
  value=f'{int(a.same_direction_patients)}/{int(a.n_pairs)}' if a.status=='DONE' else '—'
  ax.text(x,y,value,fontsize=10.8,color=INK if a.status=='DONE' else '#ABB9B1',ha='center',va='center')
legend=f.add_axes([.14,.182,.80,.028]);legend.set_xlim(0,12);legend.set_ylim(0,1);legend.axis('off')
for x,c,l in [(0,UP,'升高'),(1.5,DOWN,'降低')]:legend.scatter(x+.1,.5,s=75,c=c);legend.text(x+.32,.5,l,fontsize=10.5,color=M,va='center')
legend.scatter(3.05,.5,s=75,c=M);legend.text(3.28,.5,'实心 P < 0.05',fontsize=10.5,color=M,va='center')
legend.scatter(5.7,.5,s=75,facecolors='white',edgecolors=M,lw=1.3);legend.text(5.93,.5,'空心 P ≥ 0.05',fontsize=10.5,color=M,va='center')
legend.text(8.65,.5,'—',fontsize=15,color='#ABB9B1',va='center');legend.text(9.03,.5,'候选池未覆盖',fontsize=10.5,color=M,va='center')
line(.05,.95,.15)
text(.05,.116,'同向患者数按本行 ↑ / ↓ 统计，与各队列是否显著无关。',12,G,True)
text(.05,.078,'GBM 为非配对比较，不计算患者配对一致性；ccRCC 先汇总患者内多区域，各队列人数不合并。',10.5,M)
text(.05,.041,'范围：当前直接关系候选池，非全转录组；分母含升、降及不变的可评估配对。原排名、效应与 P 值不变。',9.2,M)
text(.95,.041,'08',14,G,True,ha='right')
f.canvas.draw();rend=f.canvas.get_renderer()
for t in f.texts:
 b=t.get_window_extent(rend);assert b.x0>=0 and b.y0>=0 and b.x1<=f.bbox.width and b.y1<=f.bbox.height
name='08_候选基因同向与患者一致性'
for ext in ['png','pdf','svg']:
 p=OUT/f'{name}.{ext}';f.savefig(p,dpi=200)
 if ext=='svg':p.write_text('\n'.join(s.rstrip() for s in p.read_text(encoding='utf8').splitlines())+'\n',encoding='utf8')
 shutil.copy2(p,VIEW/p.name)
plt.close(f)
spec=dict(version='page8_patient_consistency_v3',selection='Unchanged page8 six-cancer v2 all six genes tied at four same-significant cancers',same_direction='Each gene row dominant_direction; up counts tumor-normal>0, down counts tumor-normal<0; no within-cohort direction switch',denominator='All current-pool evaluable paired patients, including zero differences; no new sample filtering',GBM='NOT_APPLICABLE_UNPAIRED; no patient pair count',ccRCC='Use original patient-region-mean paired RNA; two cohorts kept separate',source_fields=fields,other_source_fields=['effect','p_value','n','n_up','n_down','n_equal'],RNA_family='RNA_PAIRED_CURRENT for PRAD and ccRCC; CAPT sensitivity and history excluded',new_tests=0,software={'pandas':pd.__version__,'matplotlib':matplotlib.__version__})
(OUT/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf8')
(OUT/'validation.json').write_text(json.dumps(dict(status='DONE',source_hashes_checked=True,counts_sum_to_n=True,all_effects_p_and_n_match_previous=True,selected_genes=sel.gene.tolist(),display_cells=len(counts),paired_count_cells=int(counts.status.eq('DONE').sum()),GBM_never_paired=True,visual_review='PENDING'),indent=2),encoding='utf8')
(OUT/'README_CN.md').write_text('''# 第8页：候选基因同向变化与患者一致性

本轮问题：用户要求基因图右侧也加入患者一致性，与第7页代谢物图统一样式。

输入与范围：保留六癌七队列、ABHD6/ASNS/SHMT2/SOAT1/LYPLA1/SLC7A11六个并列基因、原排序、全部点的方向及显著性。BRCA/COAD/PDAC读取原整合表配对方向计数；PRAD与ccRCC从固定Git版本的RNA_PAIRED_CURRENT表补充，不混入CAPT或历史检验族。逐项核对原效应、P、n及升降不变人数之和。

实际结果：patient_direction_counts.tsv记录每个基因每个队列的升降不变人数，以及与行主方向一致的患者数。LYPLA1为上调主方向，依次为BRCA34/45、COAD26/33、PDAC8/11、PRAD33/43、ccRCC3 5/17、ccRCC4 3/12。右侧统计不是逐患者显著检验，不能解释为这些患者各自达到P<0.05。

新手解释：同向按本行箭头计算，而非改成各癌自己的方向。ABHD6主方向向下，其计数为患者中降低的人数；其余五基因主方向向上，其计数为升高人数。反向癌种仍数与本行主方向一致的患者。因此LYPLA1的ccRCC计数较低，如实保留反证。分母包含原结果中所有可评估配对（包括不变者），不因P值或方向筛选患者。

限制/反证：候选池未覆盖不等于基因未表达或无差异。GBM非配对，没有同类患者一致性，不增设其人数列，图下注明。ccRCC先沿用患者内区域均值，不把区域数当患者数、也不跨队列合并人数。各癌种采用原RNA处理尺度和样本资格；只复用公开汇总，无原始矩阵导出、无新检验。

当前决定：PNG/PDF/SVG交付，不生成PPT。LYPLA1在ccRCC的反向及PDAC未显著均保留。

下一步：结合第9页LYPLA1四癌患者配对图讲清癌种方向与患者一致性。

复现：python code/ppt_page8_patient_consistency_v3.py。来源Git版本、路径和哈希见source_manifest.tsv。
''',encoding='utf8')
for n in ['README_CN.md','patient_direction_counts.tsv','display_genes.tsv']:shutil.copy2(OUT/n,VIEW/n)
(VIEW/'index.html').write_text(f'<!doctype html><meta charset="utf-8"><title>CAMP 第8页 · 患者一致性</title><style>body{{background:#e8eeea;margin:0;font-family:system-ui}}main{{max-width:1600px;margin:24px auto}}img{{width:100%;display:block;box-shadow:0 8px 25px #163e3320}}a{{color:#00553b}}</style><main><h2>第8页：候选基因同向变化与患者一致性</h2><img src="{name}.png"><p><a href="{name}.pdf">矢量 PDF</a> · <a href="patient_direction_counts.tsv">各队列患者计数</a> · <a href="README_CN.md">口径说明</a></p></main>',encoding='utf8')
print(counts.pivot(index='gene',columns='cohort',values='same_direction_patients').to_string());print(VIEW)
