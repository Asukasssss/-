"""Presentation-only redesign. All source rows reused unchanged from pages67_v1."""
from pathlib import Path
import shutil,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle,FancyBboxPatch

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'results/BRCA/07_INTEGRATION/20260929T120000Z_pages67_v1'
OUT=ROOT/'results/BRCA/07_INTEGRATION/20260929T140000Z_pages67_design_v2'
VIEW=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_pages_06_07_v2')
OUT.mkdir(parents=True,exist_ok=True);VIEW.mkdir(parents=True,exist_ok=True)
cov=pd.read_csv(SRC/'06_overview_data.tsv',sep='\t')
sel=pd.read_csv(SRC/'07_selected_molecules.tsv',sep='\t')
cells=pd.read_csv(SRC/'07_selected_cells.tsv',sep='\t')
for n in ['06_overview_data.tsv','07_selected_molecules.tsv','07_selected_cells.tsv','source_manifest.tsv']:
 shutil.copy2(SRC/n,OUT/n)
for fp in ['C:/Windows/Fonts/msyh.ttc','C:/Windows/Fonts/msyhbd.ttc']:font_manager.fontManager.addfont(fp)
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
G='#00553B';INK='#183D34';M='#6A7973';UP='#C77462';DOWN='#468B9E';LINE='#DDE7E0';PALE='#F3F7F3'
logo=OUT/'sysu_logo.png'
if not logo.exists():shutil.copy2(Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_SYSU_style/image2.png'),logo)
def text(f,x,y,s,size=15,col=INK,bold=False,ha='left',va='center'):
 return f.text(x,y,s,fontsize=size,color=col,weight='bold' if bold else 'normal',ha=ha,va=va)
def line(f,x1,x2,y,col=LINE,lw=.8):f.add_artist(plt.Line2D([x1,x2],[y,y],transform=f.transFigure,color=col,lw=lw))
def patch(f,x,y,w,h,c):f.add_artist(Rectangle((x,y),w,h,fc=c,ec='none',transform=f.transFigure,zorder=-5))
def base(no,title,subtitle):
 f=plt.figure(figsize=(16,9),facecolor='white')
 ax=f.add_axes([.044,.874,.177,.106]);ax.imshow(plt.imread(logo));ax.axis('off')
 text(f,.95,.928,'CAMP  /  组会汇报',10,M,ha='right');line(f,.05,.95,.876);line(f,.05,.103,.876,G,2.2)
 text(f,.05,.813,title,31,G,True);text(f,.052,.755,subtitle,13.5,M)
 text(f,.95,.038,f'{no:02}',14,G,True,ha='right')
 return f
def save(f,name):
 # Check every figure-level text rectangle fits the canvas.
 f.canvas.draw();rend=f.canvas.get_renderer();bad=[]
 for t in f.texts:
  b=t.get_window_extent(rend)
  if b.x0<0 or b.y0<0 or b.x1>f.bbox.width or b.y1>f.bbox.height:bad.append(t.get_text())
 assert not bad,bad
 for ext in ['png','pdf','svg']:
  p=OUT/f'{name}.{ext}';f.savefig(p,dpi=200)
  if ext=='svg':p.write_text('\n'.join(s.rstrip() for s in p.read_text(encoding='utf8').splitlines())+'\n',encoding='utf8')
  shutil.copy2(p,VIEW/p.name)
 plt.close(f)

f=base(6,'从代谢异常，走向候选基因','六癌七队列概况  ·  差异代谢特征 → 已映射特征 → 候选基因')
# Main plot: diverging bars make even the smallest segments readable.
text(f,.05,.691,'01',12,G,True);text(f,.085,.691,'代谢变化',18,G,True)
text(f,.720,.691,'02',12,G,True);text(f,.754,.691,'生化映射',18,G,True)
line(f,.05,.666,.661);line(f,.72,.95,.661)
text(f,.050,.625,'癌种 / 可评估特征',10.5,M)
text(f,.295,.625,'降低',12,DOWN,True,ha='center');text(f,.51,.625,'升高',12,UP,True,ha='center')
text(f,.762,.625,'已映射特征',11,M,ha='center');text(f,.906,.625,'候选基因',11,M,ha='center')
ys=np.linspace(.565,.222,7)
zero=.411;scale=.225/300
for i,r in enumerate(cov.itertuples()):
 y=ys[i]
 if i%2==0:patch(f,.047,y-.025,.91,.05,'#F7F9F7')
 text(f,.055,y+.006,r.cohort,16,G,True)
 text(f,.055,y-.016,f'{r.evaluable} 项可评估',9.5,M)
 patch(f,zero-r.down*scale,y-.010,r.down*scale,.022,DOWN)
 patch(f,zero,y-.010,r.up*scale,.022,UP)
 text(f,zero-r.down*scale-.009,y+.001,str(r.down),13,DOWN,True,ha='right')
 text(f,zero+r.up*scale+.009,y+.001,str(r.up),13,UP,True)
 text(f,.762,y+.003,str(r.mapped_features),19,INK,True,ha='center')
 text(f,.906,y+.003,str(r.candidate_genes),19,G,True,ha='center')
 text(f,.834,y+.003,'→',15,'#B2C5B7',ha='center')
f.add_artist(plt.Line2D([zero,zero],[.199,.592],transform=f.transFigure,color='#CDDAD2',lw=.9))
for n in [-300,-150,0,150,300]:
 x=zero+n*scale;text(f,x,.176,str(abs(n)),10,M,ha='center')
text(f,zero,.144,'P < 0.05 的代谢特征数',11,M,ha='center')
text(f,.720,.171,'仅统计有直接生化关系的候选',10.5,M)
line(f,.05,.95,.112)
text(f,.05,.076,'队列覆盖与映射范围不同，数量不代表癌种异常程度。',12,INK)
text(f,.05,.037,'名义 P < 0.05；GBM 非配对，其余按患者配对；ccRCC 两队列分列。',9.5,M)
save(f,'06_代谢异常到候选基因_v2')

f=base(7,'跨癌代谢变化：共性与方向差异','10 个代表代谢物  ·  同一化学标识逐项比较，保留不显著与未覆盖结果')
labels={'KEGG:C00328':'犬尿氨酸','KEGG:C00137':'肌醇','KEGG:C02990':'棕榈酰肉碱','KEGG:C00319':'鞘氨醇','KEGG:C00025':'谷氨酸','KEGG:C00031':'葡萄糖','KEGG:C00065':'丝氨酸','KEGG:C00123':'亮氨酸','KEGG:C00864':'泛酸（维生素 B5）','KEGG:C00047':'赖氨酸'}
order=cov.cohort.tolist()
# Scientific point matrix. Plot coordinates are shared by points and labels.
ax=f.add_axes([.05,.198,.90,.505]);ax.set_xlim(-3.8,7.5);ax.set_ylim(11.25,-1.05);ax.axis('off')
xx=np.arange(7)+.55
yy=[.35,1.25,2.15,3.05,5.1,6,6.9,7.8,8.7,9.6]
for j,c in enumerate(order):ax.text(xx[j],-.53,c,ha='center',va='center',fontsize=12.5,color=G,weight='bold')
ax.axhline(-.08,color=LINE,lw=.9)
for j in range(7):ax.plot([xx[j],xx[j]],[.02,10.13],c='#EAF0EB',lw=.7,zorder=-3)
# Group labels are horizontal and distinct from row names.
ax.text(-3.7,.35,'同向\n趋势',ha='left',va='center',fontsize=13,color=G,weight='bold',linespacing=1.5)
ax.text(-3.7,5.1,'相反\n方向',ha='left',va='center',fontsize=13,color=G,weight='bold',linespacing=1.5)
ax.plot([-3.7,7.45],[4.08,4.08],c=LINE,lw=1)
for i,r in enumerate(sel.itertuples()):
 y=yy[i]
 if i in [0,5]:ax.add_patch(Rectangle((-2.8,y-.40),10.25,.80,fc=PALE,ec='none',zorder=-4))
 ax.text(-2.72,y,labels[r.key],va='center',ha='left',fontsize=13.5,color=INK,weight='bold' if i in [0,5] else 'normal')
 for j,c in enumerate(order):
  d=cells[(cells.key==r.key)&(cells.cohort==c)].iloc[0]
  if d.status=='EVALUABLE':
   col=UP if d.effect>0 else DOWN;ax.scatter(xx[j],y,s=170,facecolors=col if d.p<.05 else 'white',edgecolors=col,linewidths=1.55,zorder=5)
  else:ax.text(xx[j],y,{'ABSENT_IN_SOURCE':'—','NOT_EVALUABLE':'×','AMBIGUOUS_MULTIPLE_FEATURES':'?'}[d.status],fontsize=17,ha='center',va='center',color='#B1BDB6')
# Compact legend on one aligned line, below the graph.
legend=f.add_axes([.24,.155,.68,.033]);legend.set_xlim(0,10);legend.set_ylim(0,1);legend.axis('off')
for x,col,label in [(0,UP,'升高'),(1.6,DOWN,'降低')]:legend.scatter(x+.1,.5,s=90,c=col);legend.text(x+.32,.5,label,va='center',fontsize=10.5,color=M)
legend.scatter(3.3,.5,s=90,c=M);legend.text(3.53,.5,'实心 P < 0.05',va='center',fontsize=10.5,color=M)
legend.scatter(5.7,.5,s=90,facecolors='white',edgecolors=M,linewidths=1.3);legend.text(5.93,.5,'空心 P ≥ 0.05',va='center',fontsize=10.5,color=M)
legend.text(7.8,.5,'—',fontsize=15,color='#B1BDB6',va='center');legend.text(8.15,.5,'未覆盖',fontsize=10.5,color=M,va='center')
line(f,.05,.95,.132)
text(f,.05,.103,'犬尿氨酸',12.5,G,True);text(f,.155,.103,'6 个可评估队列均升高，4 个显著',11.5,INK)
text(f,.56,.103,'葡萄糖',12.5,G,True);text(f,.635,.103,'COAD / PRAD 降低，GBM / ccRCC 升高',11.5,INK)
text(f,.05,.066,'同向不等于各队列均显著；相反方向提供线索，尚未检验癌种间差异。',10.5,M)
text(f,.05,.032,'代表分子沿用同一覆盖度排序规则；圆点等大，颜色不表示效应幅度。完整数据与选图规则见附表。',9.2,M)
save(f,'07_跨癌代谢方向点阵_v2')

# Verify no quantitative input changed in this layout-only revision.
validation={'data_unchanged':True,'no_new_tests':True,'visualization_changes':['diverging bars with every count visible','two-column mapping coverage','equal-size filled/hollow dot matrix','original template logo','reduced prose and larger text'],'no_ppt_created':True}
for n in ['06_overview_data.tsv','07_selected_molecules.tsv','07_selected_cells.tsv']:
 assert (OUT/n).read_bytes()==(SRC/n).read_bytes();validation[n]=hashlib.sha256((OUT/n).read_bytes()).hexdigest()
(OUT/'validation.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2),encoding='utf8')
(OUT/'README_CN.md').write_text('''# 第6–7页版式重设计 v2

只调整绘图与排版，输入逐字节复用v1。没有重新选择分子、改变P值或筛选范围，没有生成PPT。

第6页：升降分列，以0为中心的条形图，每个数值均在条形外标出；非显著背景通过可评估分母保留。直接关系数移至附表，正文聚焦已映射分子与候选基因。

第7页：原10个分子和70个格子保留。等大圆点，红升蓝降；实心P<0.05，空心P>=0.05；破折号为当前表未覆盖。本次选中格子没有多特征身份歧义或不可评估状态，因此只展示实际用到的图例。同向趋势不表示每位患者同向或所有队列均显著。

范围：GBM非配对，BRCA45、COAD33、PDAC11、PRAD43、ccRCC3 17、ccRCC4 12对；原始方法和选图规则见v1分析规格。癌种分化为描述性方向对照，不是已完成癌种交互检验。

校徽来自用户提供的Pre_IBD模板提取资产；科研图表均由代码绘制。复现：python code/ppt_pages67_v2.py。
''',encoding='utf8')
for n in ['README_CN.md','validation.json','06_overview_data.tsv','07_selected_cells.tsv']:shutil.copy2(OUT/n,VIEW/n)
html='<!doctype html><meta charset="utf-8"><title>CAMP 第6–7页 · v2</title><style>body{background:#e8eeea;margin:0;font-family:system-ui}main{max-width:1440px;margin:24px auto}img{width:100%;display:block;margin:25px 0;box-shadow:0 8px 25px #163e3320}a{color:#00553b}</style><main><h2>CAMP · 第6–7页重设计</h2>'
for n in ['06_代谢异常到候选基因_v2','07_跨癌代谢方向点阵_v2']:html+=f'<img src="{n}.png"><a href="{n}.pdf">矢量 PDF</a>'
html+='<p>数据与上一版一致，仅调整版式。</p></main>';(VIEW/'index.html').write_text(html,encoding='utf8')
print(VIEW)
