"""Page 7: rank same-direction nominal significance across cancer types."""
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
SRC=ROOT/'results/BRCA/07_INTEGRATION/20260929T120000Z_pages67_v1'
PREV=ROOT/'results/BRCA/07_INTEGRATION/20260929T140000Z_pages67_design_v2'
OUT=ROOT/'results/BRCA/07_INTEGRATION/20260929T150000Z_page7_same_direction_v3'
VIEW=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_page07_same_direction_v3')
OUT.mkdir(parents=True,exist_ok=True);VIEW.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(SRC/'07_all_identity_matrix.tsv',sep='\t');rows=[]
for k,z in d.groupby('key'):
 if k in ['KEGG:C00031','KEGG:C00157']:continue # user excludes glucose; generic phosphatidylcholine is not one molecule
 z=z[z.status=='EVALUABLE'].copy();z['cancer']=z.cohort.replace({'ccRCC3':'ccRCC','ccRCC4':'ccRCC'})
 s=z[z.p<.05];up=s[s.effect>0];dn=s[s.effect<0];u=up.cancer.nunique();v=dn.cancer.nunique()
 direction='up' if u>v else 'down' if v>u else 'tie';main=up if direction=='up' else dn
 rows.append(dict(key=k,source_names=' | '.join(z.source_names.dropna().unique()),same_significant_cancers=max(u,v),opposite_significant_cancers=min(u,v),dominant_direction=direction,evaluable_cancers=z.cancer.nunique(),same_significant_cohorts=len(main),evaluable_cohorts=len(z)))
ranking=pd.DataFrame(rows).sort_values(['same_significant_cancers','opposite_significant_cancers','evaluable_cancers','same_significant_cohorts','key'],ascending=[False,True,False,False,True])
sel=ranking[(ranking.dominant_direction!='tie')&(ranking.same_significant_cancers>=3)].head(8).copy()
assert len(sel)==8 and sel.iloc[0]['key']=='KEGG:C00791'
display=d[d.key.isin(sel.key)].copy();assert len(display)==56
ranking.to_csv(OUT/'ranking_all.tsv',sep='\t',index=False);sel.to_csv(OUT/'display_molecules.tsv',sep='\t',index=False);display.to_csv(OUT/'display_cells.tsv',sep='\t',index=False)
shutil.copy2(SRC/'source_manifest.tsv',OUT/'source_manifest.tsv');shutil.copy2(PREV/'sysu_logo.png',OUT/'sysu_logo.png')
for fp in ['C:/Windows/Fonts/msyh.ttc','C:/Windows/Fonts/msyhbd.ttc']:font_manager.fontManager.addfont(fp)
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
G='#00553B';INK='#183D34';M='#6A7973';UP='#C77462';DOWN='#468B9E';LINE='#DDE7E0';PALE='#F3F7F3'
f=plt.figure(figsize=(16,9),facecolor='white')
def text(x,y,s,size=15,col=INK,bold=False,ha='left'):return f.text(x,y,s,fontsize=size,color=col,weight='bold' if bold else 'normal',ha=ha,va='center')
def line(x1,x2,y,c=LINE,w=.8):f.add_artist(plt.Line2D([x1,x2],[y,y],transform=f.transFigure,color=c,lw=w))
ax=f.add_axes([.044,.874,.177,.106]);ax.imshow(plt.imread(OUT/'sysu_logo.png'));ax.axis('off')
text(.95,.928,'CAMP  /  组会汇报',10,M,ha='right');line(.05,.95,.876);line(.05,.103,.876,G,2.2)
text(.05,.813,'哪些代谢物在更多癌种中同向改变？',30,G,True)
text(.052,.755,'按同一方向 P < 0.05 的癌种数排序  ·  ccRCC 两队列只计一个癌种',13.5,M)
labels={'KEGG:C00791':'肌酐','KEGG:C02630':'2-羟基戊二酸¹','KEGG:C02862':'丁酰肉碱','KEGG:C00294':'肌苷','KEGG:C00037':'甘氨酸','KEGG:C00208':'麦芽糖','KEGG:C00328':'犬尿氨酸','KEGG:C02712':'N-乙酰甲硫氨酸'}
assert set(sel.key)==set(labels)
order=['BRCA','COAD','GBM','PDAC','PRAD','ccRCC3','ccRCC4']
ax=f.add_axes([.05,.218,.90,.493]);ax.set_xlim(-3.3,9.05);ax.set_ylim(8.2,-1.2);ax.axis('off')
xx=np.arange(7)+.1
for j,c in enumerate(order):ax.text(xx[j],-.58,c,ha='center',va='center',fontsize=12,color=G,weight='bold')
ax.text(8.05,-.66,'同向显著 / 可评估',fontsize=10.7,ha='center',va='center',color=G,weight='bold');ax.text(8.05,-.29,'癌种数',fontsize=10,ha='center',va='center',color=M)
ax.plot([-3.3,9.05],[-.07,-.07],color=LINE,lw=.8)
for j in range(7):ax.plot([xx[j],xx[j]],[.02,7.8],c='#ECF0EC',lw=.6,zorder=-3)
for i,r in enumerate(sel.itertuples()):
 y=i+.35
 if i%2==0:ax.add_patch(Rectangle((-3.3,y-.43),12.35,.86,fc='#F7F9F7',ec='none',zorder=-4))
 ax.text(-3.2,y,labels[r.key],fontsize=13.5,color=INK,va='center',weight='bold' if i<4 else 'normal')
 for j,c in enumerate(order):
  a=display[(display.key==r.key)&(display.cohort==c)].iloc[0]
  if a.status=='EVALUABLE':
   color=UP if a.effect>0 else DOWN;ax.scatter(xx[j],y,s=155,facecolors=color if a.p<.05 else 'white',edgecolors=color,linewidths=1.5,zorder=3)
  else:ax.text(xx[j],y,{'ABSENT_IN_SOURCE':'—','AMBIGUOUS_MULTIPLE_FEATURES':'?','NOT_EVALUABLE':'×'}[a.status],ha='center',va='center',color='#B1BDB6',fontsize=17)
 color=UP if r.dominant_direction=='up' else DOWN
 ax.text(7.4,y,'↑' if r.dominant_direction=='up' else '↓',fontsize=18,color=color,ha='center',va='center')
 ax.text(8.05,y,f'{r.same_significant_cancers} / {r.evaluable_cancers}',fontsize=19,color=color,ha='center',va='center',weight='bold')
ax.plot([6.9,6.9],[-.94,7.8],color=LINE,lw=.8)
legend=f.add_axes([.19,.173,.75,.031]);legend.set_xlim(0,11);legend.set_ylim(0,1);legend.axis('off')
for x,c,l in [(0,UP,'升高'),(1.6,DOWN,'降低')]:legend.scatter(x+.1,.5,s=85,c=c);legend.text(x+.32,.5,l,fontsize=10.5,color=M,va='center')
legend.scatter(3.3,.5,s=85,c=M);legend.text(3.53,.5,'实心 P < 0.05',fontsize=10.5,color=M,va='center')
legend.scatter(6,.5,s=85,facecolors='white',edgecolors=M,linewidths=1.3);legend.text(6.23,.5,'空心 P ≥ 0.05',fontsize=10.5,color=M,va='center')
legend.text(8.9,.5,'—',fontsize=15,color='#B1BDB6',va='center');legend.text(9.3,.5,'未覆盖',fontsize=10.5,color=M,va='center')
line(.05,.95,.15)
text(.05,.116,'肌酐：5 癌显著降低',14,G,True);text(.42,.116,'同向次数最多，不要求所有队列方向一致。',12,INK)
text(.05,.078,'肌苷在 ccRCC3 显著升高，反向结果保留；未显著方向不计入同向显著癌种数。',10.5,M)
text(.05,.041,'¹ 2-羟基戊二酸沿用源注释，未区分 D/L。并列排序及完整来源见附表；GBM 为非配对比较。',9.3,M)
text(.95,.041,'07',14,G,True,ha='right')
f.canvas.draw();rend=f.canvas.get_renderer()
for t in f.texts:
 b=t.get_window_extent(rend);assert b.x0>=0 and b.y0>=0 and b.x1<=f.bbox.width and b.y1<=f.bbox.height
name='07_跨癌同向重复最多_v3'
for ext in ['png','pdf','svg']:
 p=OUT/f'{name}.{ext}';f.savefig(p,dpi=200)
 if ext=='svg':p.write_text('\n'.join(s.rstrip() for s in p.read_text(encoding='utf8').splitlines())+'\n',encoding='utf8')
 shutil.copy2(p,VIEW/p.name)
plt.close(f)
spec={'selection':'Top 8 by same-direction nominally significant cancer count; ties by fewer reverse significant cancers, more evaluable cancers, more same-direction significant cohorts, chemical key','p':.05,'ccRCC_count':'Either cohort significant in a direction counts ccRCC once for that direction; if significant opposite directions occur, record ccRCC in both direction-specific counts; never sum the two direction-specific counts as distinct cancer total','excluded':['KEGG:C00031 user requested removal','KEGG:C00157 generic lipid class'],'no_new_patient_tests':True,'identity':'Inherited KEGG/HMDB source keys, unique evaluable feature per cohort; no alias or peak pooling; no new chemical identification','source_matrix_sha256':hashlib.sha256((SRC/'07_all_identity_matrix.tsv').read_bytes()).hexdigest(),'display_count':8,'cell_count':56}
(OUT/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf8')
(OUT/'README_CN.md').write_text('''# 第7页：跨癌同向显著重复 v3

用户要求仅展示同向最多的几个。本版删除先分同向/反向两组的选择限制，从全部已标准化KEGG/HMDB结果中按同向P<0.05癌种数排序。最多方向为主方向；并列方向者不入主展示。先按主方向显著癌种数降序，再按反向显著癌种数升序、可评估癌种数降序、同向显著队列数降序、化学键字典序。展示前8个，葡萄糖保持排除。

ccRCC3/4分别绘图，在任一方向达到P<0.05即为该方向贡献一个ccRCC癌种计数，最多一次；如两队列各有显著相反方向，两方向计数均记录，但不得相加当不同癌种。本次展示8个分子没有这种双显著ccRCC冲突；肌苷ccRCC3显著上升、ccRCC4不显著下降，故计一个反向癌种。覆盖分母为至少一个队列可唯一对应且可评估的癌种数。

肌酐5癌下降；2-羟基戊二酸、丁酰肉碱4癌上升；肌苷4癌下降且1癌反向；甘氨酸、麦芽糖、犬尿氨酸、N-乙酰甲硫氨酸3癌上升。2-HG源注释未区分D/L，不外推构型。此页不是P值最小排名，也不要求全部队列同向。空心、反向与未覆盖均完整展示。

原P值和效应全部复用，不重算统计，不生成PPT。第6页未修改。复现：python code/ppt_page7_same_direction_v3.py。
''',encoding='utf8')
for n in ['display_molecules.tsv','display_cells.tsv','ranking_all.tsv','analysis_spec.json','README_CN.md']:shutil.copy2(OUT/n,VIEW/n)
(VIEW/'index.html').write_text(f'<!doctype html><meta charset="utf-8"><title>CAMP 第7页 · 同向重复</title><style>body{{background:#e8eeea;margin:0;font-family:system-ui}}main{{max-width:1440px;margin:24px auto}}img{{width:100%;display:block;margin:24px 0;box-shadow:0 8px 25px #163e3320}}a{{color:#00553b}}</style><main><h2>第7页：同向显著重复最多的代表</h2><img src="{name}.png"><a href="{name}.pdf">矢量 PDF</a> · <a href="display_molecules.tsv">排序与癌种数</a><p>按同一方向P&lt;0.05癌种数排序；葡萄糖已排除。</p></main>',encoding='utf8')
print(sel.to_string(index=False));print(VIEW)
