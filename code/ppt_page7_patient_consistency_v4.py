"""Add per-cohort paired patient direction counts to the frozen page7 selection."""
from pathlib import Path
import json,hashlib,shutil
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'results/BRCA/07_INTEGRATION/20260928T150000Z_ppt_redraw_v1'
PREV=ROOT/'results/BRCA/07_INTEGRATION/20260929T150000Z_page7_same_direction_v3'
OUT=ROOT/'results/BRCA/07_INTEGRATION/20260929T190000Z_page7_patient_consistency_v4'
VIEW=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_page07_patient_consistency_v4')
OUT.mkdir(parents=True,exist_ok=True);VIEW.mkdir(parents=True,exist_ok=True)
sel=pd.read_csv(PREV/'display_molecules.tsv',sep='\t');d=pd.read_csv(PREV/'display_cells.tsv',sep='\t')
manifest=pd.read_csv(OLD/'source_manifest.tsv',sep='\t');sources=[];tabs={}
for c in ['brca','coad','pdac','prad','ccrcc']:
 r=manifest[manifest.id.eq(c+'_met')].iloc[0].to_dict();p=OLD/'sources'/f'{c}_met.tsv';assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'];sources.append(r);tabs[c]=pd.read_csv(p,sep='\t')
for name in ['display_cells.tsv','display_molecules.tsv']:
 p=PREV/name;sources.append(dict(id='previous_'+name,git_path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size,git_commit='33e2083'))
records=[]
order=['BRCA','COAD','GBM','PDAC','PRAD','ccRCC3','ccRCC4'];paired=[c for c in order if c!='GBM']
for r in sel.itertuples():
 for c in order:
  a=d[d.key.eq(r.key)&d.cohort.eq(c)].iloc[0]
  base=dict(key=r.key,cohort=c,main_direction=r.dominant_direction,cohort_effect=a.effect,cohort_p=a.p,same_direction_patients=np.nan,n_pairs=np.nan,up_patients=np.nan,down_patients=np.nan,equal_patients=np.nan,same_fraction=np.nan)
  if a.status!='EVALUABLE':records.append(dict(base,status=a.status));continue
  if c=='GBM':records.append(dict(base,status='NOT_APPLICABLE_UNPAIRED'));continue
  tab=tabs['ccrcc' if c.startswith('ccRCC') else c.lower()]
  z=tab[tab.metabolite_key.eq(r.key)]
  if c.startswith('ccRCC'):z=z[z.cohort.eq(c)]
  assert len(z)==1,(r.key,c);z=z.iloc[0]
  assert np.isclose(z.effect,a.effect,rtol=0,atol=1e-12) and np.isclose(z.p_value,a.p,rtol=1e-12,atol=1e-15)
  uc,dc,ec=('pairs_higher','pairs_lower','pairs_equal') if c in ['BRCA','COAD'] else ('n_up','n_down','n_equal')
  up,down,equal,n=[int(z[k]) for k in [uc,dc,ec,'n']];assert up+down+equal==n
  same=up if r.dominant_direction=='up' else down
  base.update(same_direction_patients=same,n_pairs=n,up_patients=up,down_patients=down,equal_patients=equal,same_fraction=same/n)
  records.append(dict(base,status='DONE'))
counts=pd.DataFrame(records);assert len(counts)==56
counts.to_csv(OUT/'patient_direction_counts.tsv',sep='\t',index=False)
sel.to_csv(OUT/'display_molecules.tsv',sep='\t',index=False);d.to_csv(OUT/'display_cells.tsv',sep='\t',index=False)
pd.DataFrame(sources).to_csv(OUT/'source_manifest.tsv',sep='\t',index=False)
shutil.copy2(PREV/'sysu_logo.png',OUT/'sysu_logo.png')
for fp in ['C:/Windows/Fonts/msyh.ttc','C:/Windows/Fonts/msyhbd.ttc']:font_manager.fontManager.addfont(fp)
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
G='#00553B';INK='#183D34';M='#6A7973';UP='#C77462';DOWN='#468B9E';LINE='#DDE7E0'
labels={'KEGG:C00791':'肌酐','KEGG:C02630':'2-羟基戊二酸¹','KEGG:C02862':'丁酰肉碱','KEGG:C00294':'肌苷','KEGG:C00037':'甘氨酸','KEGG:C00208':'麦芽糖','KEGG:C00328':'犬尿氨酸','KEGG:C02712':'N-乙酰甲硫氨酸'}
f=plt.figure(figsize=(16,9),facecolor='white')
def text(x,y,s,size=15,col=INK,bold=False,ha='left'):return f.text(x,y,s,fontsize=size,color=col,weight='bold' if bold else 'normal',ha=ha,va='center')
def line(x1,x2,y,c=LINE,w=.8):f.add_artist(plt.Line2D([x1,x2],[y,y],transform=f.transFigure,color=c,lw=w))
ax=f.add_axes([.044,.874,.177,.106]);ax.imshow(plt.imread(OUT/'sysu_logo.png'));ax.axis('off')
text(.95,.928,'CAMP  /  组会汇报',10,M,ha='right');line(.05,.95,.876);line(.05,.103,.876,G,2.2)
text(.05,.813,'跨癌同向变化，在多少患者中重复？',30,G,True)
text(.052,.755,'沿用同向显著癌种数排序  ·  右侧统计与每行主方向一致的患者数 / 可评估配对数',13,M)
# One common canvas keeps all rows and count columns exactly aligned.
ax=f.add_axes([.05,.236,.90,.472]);ax.set_xlim(0,100);ax.set_ylim(8.15,-1.5);ax.axis('off')
xx=np.linspace(23,56.5,7);cx=63.3;px=np.linspace(72.1,97.6,6)
ax.text(0,-.45,'代谢物',fontsize=11,color=G,weight='bold',va='center')
for x,c in zip(xx,order):ax.text(x,-.45,c,fontsize=10.4,color=G,weight='bold',ha='center',va='center')
ax.text(cx,-.64,'同向显著 / 可评估',fontsize=8.6,color=G,weight='bold',ha='center',va='center');ax.text(cx,-.25,'癌种数',fontsize=9,color=M,ha='center',va='center')
ax.text(84.8,-1.17,'患者一致性',fontsize=13,color=G,weight='bold',ha='center',va='center')
for x,c in zip(px,paired):ax.text(x,-.45,c,fontsize=9.3,color=G,weight='bold',ha='center',va='center')
ax.plot([0,100],[.02,.02],c=LINE,lw=.8)
for x in [59.7,68.4]:ax.plot([x,x],[-1.4,7.8],c=LINE,lw=.8)
for i,r in enumerate(sel.itertuples()):
 y=i+.48;col=UP if r.dominant_direction=='up' else DOWN
 if i%2==0:ax.add_patch(Rectangle((0,y-.44),100,.88,fc='#F5F8F5',ec='none',zorder=-3))
 ax.text(.6,y,labels[r.key],fontsize=12.3,color=INK,weight='bold' if i<4 else 'normal',va='center')
 for x,c in zip(xx,order):
  a=d[d.key.eq(r.key)&d.cohort.eq(c)].iloc[0]
  if a.status=='EVALUABLE':
   color=UP if a.effect>0 else DOWN;ax.scatter(x,y,s=120,facecolors=color if a.p<.05 else 'white',edgecolors=color,lw=1.4)
  else:ax.text(x,y,'—',fontsize=15,color='#ABB9B1',ha='center',va='center')
 ax.text(cx-2.4,y,'↑' if r.dominant_direction=='up' else '↓',fontsize=16,color=col,ha='center',va='center')
 ax.text(cx+.5,y,f'{r.same_significant_cancers}/{r.evaluable_cancers}',fontsize=16,color=col,ha='center',va='center',weight='bold')
 for x,c in zip(px,paired):
  a=counts[counts.key.eq(r.key)&counts.cohort.eq(c)].iloc[0]
  value=f'{int(a.same_direction_patients)}/{int(a.n_pairs)}' if a.status=='DONE' else '—'
  ax.text(x,y,value,fontsize=10.8,color=INK if a.status=='DONE' else '#ABB9B1',ha='center',va='center')
legend=f.add_axes([.14,.182,.80,.028]);legend.set_xlim(0,12);legend.set_ylim(0,1);legend.axis('off')
for x,c,l in [(0,UP,'升高'),(1.5,DOWN,'降低')]:legend.scatter(x+.1,.5,s=75,c=c);legend.text(x+.32,.5,l,fontsize=10.5,color=M,va='center')
legend.scatter(3.05,.5,s=75,c=M);legend.text(3.28,.5,'实心 P < 0.05',fontsize=10.5,color=M,va='center')
legend.scatter(5.7,.5,s=75,facecolors='white',edgecolors=M,lw=1.3);legend.text(5.93,.5,'空心 P ≥ 0.05',fontsize=10.5,color=M,va='center')
legend.text(8.65,.5,'—',fontsize=15,color='#ABB9B1',va='center');legend.text(9.03,.5,'未覆盖',fontsize=10.5,color=M,va='center')
line(.05,.95,.15)
text(.05,.116,'同向患者数按本行 ↑ / ↓ 统计，与各队列是否显著无关。',12,G,True)
text(.05,.078,'GBM 为非配对比较，不计算患者配对一致性；ccRCC 先汇总患者内多区域，各队列人数不合并。',10.5,M)
text(.05,.041,'分母包含升高、降低与不变的可评估配对；使用原处理矩阵。¹ 2-羟基戊二酸未区分 D/L。原排序、效应与 P 值不变。',9.2,M)
text(.95,.041,'07',14,G,True,ha='right')
f.canvas.draw();rend=f.canvas.get_renderer()
for t in f.texts:
 b=t.get_window_extent(rend);assert b.x0>=0 and b.y0>=0 and b.x1<=f.bbox.width and b.y1<=f.bbox.height
name='07_代谢物同向与患者一致性'
for ext in ['png','pdf','svg']:
 p=OUT/f'{name}.{ext}';f.savefig(p,dpi=200)
 if ext=='svg':p.write_text('\n'.join(s.rstrip() for s in p.read_text(encoding='utf8').splitlines())+'\n',encoding='utf8')
 shutil.copy2(p,VIEW/p.name)
plt.close(f)
spec=dict(version='page7_patient_consistency_v4',selection='Unchanged page7 v3 eight molecules; no reranking',same_direction='Each selected row dominant_direction: positive tumor-normal for up, negative for down; includes P>=0.05 and opposing cohorts',denominator='All evaluable paired patients including zero differences on unchanged author processed matrices; not restricted to observed before imputation',cohort_pooling='None; ccRCC cohorts displayed separately; multiple regions already averaged per patient in original results',GBM='NOT_APPLICABLE_UNPAIRED; no artificial pairing or comparison to reference mean',new_statistics='Only descriptive count column reuse; no inferential tests, matrix access, normalization or P changes',source_count_columns={'BRCA/COAD':['pairs_higher','pairs_lower','pairs_equal'],'PDAC/PRAD/ccRCC':['n_up','n_down','n_equal']},software={'pandas':pd.__version__,'matplotlib':matplotlib.__version__})
(OUT/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf8')
(OUT/'validation.json').write_text(json.dumps(dict(status='DONE',source_hashes_checked=True,all_counts_sum_to_n=True,all_effects_and_p_match_previous=True,display_molecules=8,display_cells=56,paired_counts=int(counts.status.eq('DONE').sum()),GBM_never_paired=True,visual_review='PENDING'),indent=2),encoding='utf8')
(OUT/'README_CN.md').write_text('''# 第7页新增患者一致性

本轮问题：用户要求右侧增加多少患者同向。保留上一版8个代谢物、原有7队列效应方向及显著性、同向显著癌种计数，右侧新增6个有配对队列的患者计数。

输入与范围：读取既有公开汇总表，BRCA/COAD使用pairs_higher/lower/equal；PDAC/PRAD/ccRCC使用n_up/down/equal。逐行核对化学键唯一、队列、原效应、P、n及计数之和。没有重新访问患者矩阵或重算统计。来源哈希与Git版本见source_manifest.tsv。

实际结果：patient_direction_counts.tsv保留每分子每队列的升/降/不变、同向人数及分母。肌酐按降低统计为BRCA37/45、COAD32/33、PDAC5/11、PRAD30/43、ccRCC3 16/17、ccRCC4 11/12。

新手解释：同向是与本行箭头的主方向相同，不是与每个队列自己的方向相同。例如肌酐主方向降低，PDAC虽平均升高，仍数其中降低的5/11患者。显著与不显著队列都展示；患者计数不做P值过滤。右侧不是“每位患者均统计显著”的人数。分母为该特征可评估的全部配对，包括差值为0者；零差值不计入同向分子。

限制/反证：使用原作者处理矩阵，部分值可能有原作者填补，不等于全部原始实测检出。GBM无患者配对，不能计算同一种配对一致性，右侧不设GBM计数列并在图下注明。ccRCC先在患者内汇总区域，不能把区域当患者；两个队列人数也不合并。未覆盖保持横线，不当0。原来的显著反向和不显著方向不删。

当前决定：图片、PDF和SVG交付，不生成PPT；癌种同向数与患者同向数分开呈现。

下一步：以统一口径解读癌种和患者两层重复性。

复现：python code/ppt_page7_patient_consistency_v4.py。
''',encoding='utf8')
for n in ['README_CN.md','patient_direction_counts.tsv','display_molecules.tsv']:shutil.copy2(OUT/n,VIEW/n)
(VIEW/'index.html').write_text(f'<!doctype html><meta charset="utf-8"><title>CAMP 第7页 · 患者一致性</title><style>body{{background:#e8eeea;margin:0;font-family:system-ui}}main{{max-width:1600px;margin:24px auto}}img{{width:100%;display:block;box-shadow:0 8px 25px #163e3320}}a{{color:#00553b}}</style><main><h2>第7页：代谢物同向变化与患者一致性</h2><img src="{name}.png"><p><a href="{name}.pdf">矢量 PDF</a> · <a href="patient_direction_counts.tsv">各队列患者计数</a> · <a href="README_CN.md">口径说明</a></p></main>',encoding='utf8')
print(counts.pivot(index='key',columns='cohort',values='same_direction_patients').to_string());print(VIEW)
