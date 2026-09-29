"""Motivation slide: why follow LYPLA1, using frozen four-cancer aggregates."""
from pathlib import Path
import json,hashlib,shutil
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch,Rectangle
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'results/BRCA/07_INTEGRATION/20260929T180100Z_page9_lypla1_paired_v1'
GEN=ROOT/'results/BRCA/07_INTEGRATION/20260929T200000Z_page8_patient_consistency_v3'
OUT=ROOT/'results/BRCA/07_INTEGRATION/20260929T210000Z_focus_reason_v1'
VIEW=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_LYPLA1_focus_reason_v1')
OUT.mkdir(parents=True,exist_ok=True);VIEW.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(SRC/'paired_summary.tsv',sep='\t');g=pd.read_csv(GEN/'patient_direction_counts.tsv',sep='\t')
assert d.cancer.tolist()==['BRCA','COAD','PDAC','PRAD'] and d.mean_tumor_minus_normal.gt(0).all() and d.p_value.lt(.05).sum()==3
for r in d.itertuples():
 q=g[g.gene.eq('LYPLA1')&g.cohort.eq(r.cancer)].iloc[0]
 assert q.same_direction_patients==r.up_pairs and q.n_pairs==r.n_pairs and np.isclose(q.cohort_p,r.p_value)
d.to_csv(OUT/'four_cancer_evidence.tsv',sep='\t',index=False)
inputs=[SRC/'paired_summary.tsv',SRC/'source_manifest_frozen.tsv',GEN/'patient_direction_counts.tsv',GEN/'sysu_logo.png']
pd.DataFrame([dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in inputs]).to_csv(OUT/'source_manifest.tsv',sep='\t',index=False)
shutil.copy2(SRC/'source_manifest_frozen.tsv',OUT/'source_manifest_frozen.tsv')
for fp in ['C:/Windows/Fonts/msyh.ttc','C:/Windows/Fonts/msyhbd.ttc']:font_manager.fontManager.addfont(fp)
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
G='#00553B';INK='#183D34';M='#6A7973';UP='#C77462';LINE='#DDE7E0';PALE='#F2F7F3'
f=plt.figure(figsize=(16,9),facecolor='white')
def text(x,y,s,size=15,col=INK,bold=False,ha='left'):return f.text(x,y,s,fontsize=size,color=col,weight='bold' if bold else 'normal',ha=ha,va='center')
def line(x1,x2,y,c=LINE,w=.8):f.add_artist(plt.Line2D([x1,x2],[y,y],transform=f.transFigure,color=c,lw=w))
def box(x,y,w,h,fc=PALE):f.add_artist(FancyBboxPatch((x,y),w,h,transform=f.transFigure,boxstyle='round,pad=0,rounding_size=.012',facecolor=fc,edgecolor='none',zorder=-3))
ax=f.add_axes([.044,.874,.177,.106]);ax.imshow(plt.imread(GEN/'sysu_logo.png'));ax.axis('off')
text(.95,.928,'CAMP  /  组会汇报',10,M,ha='right');line(.05,.95,.876);line(.05,.103,.876,G,2.2)
text(.05,.813,'从跨癌候选中，为什么关注 LYPLA1？',30,G,True)
text(.052,.755,'代谢关系提供候选来源，跨癌表达与患者同向性支持进一步追踪',14,M)
text(.052,.682,'发现路径',16,G,True)
cards=[(.518,'01','CAMP 代谢异常','肿瘤—参照比较，定位变化的代谢物'),(.366,'02','直接生化关系映射','从代谢物连接到相关酶、转运体基因'),(.214,'03','关注 LYPLA1','在候选池中发现跨癌重复的表达线索')]
for y,num,title,desc in cards:
 box(.05,y,.345,.117)
 text(.068,y+.079,num,17,G,True)
 text(.105,y+.078,title,17,G,True)
 text(.105,y+.036,desc,11.1,M)
 if num!='03':
  f.add_artist(plt.Line2D([.222,.222],[y-.007,y-.026],transform=f.transFigure,color='#9AB2A3',lw=1.3))
  f.text(.222,y-.030,'▼',fontsize=9,color='#9AB2A3',ha='center',va='center')
text(.447,.682,'四癌中观察到的表达线索',16,G,True)
text(.949,.683,'组织 RNA · 患者配对',10.5,M,ha='right')
heads=[(.472,'癌种'),(.590,'平均方向'),(.712,'P 值'),(.858,'上调患者 / 配对数')]
for x,s in heads:text(x,.619,s,10.8,M,ha='center')
line(.444,.95,.591)
for i,r in enumerate(d.itertuples()):
 y=.551-i*.081
 if i%2==0:box(.443,y-.035,.507,.070,'#F7F9F7')
 text(.472,y,r.cancer,14,G,True,ha='center')
 text(.590,y,'↑ 上调',13,UP,True,ha='center')
 text(.712,y,f'{r.p_value:.3g}'+(' *' if r.p_value<.05 else ''),12.5,INK,ha='center')
 text(.823,y,f'{r.up_pairs} / {r.n_pairs}',13,INK,True,ha='center')
 text(.927,y,f'{r.up_fraction:.1%}',11.5,M,ha='right')
 # Compact bar encodes the fraction; absolute counts remain explicit above.
 f.add_artist(Rectangle((.788,y-.025),.140,.004,transform=f.transFigure,fc='#E2EAE4',ec='none'))
 f.add_artist(Rectangle((.788,y-.025),.140*r.up_fraction,.004,transform=f.transFigure,fc=G,ec='none'))
text(.447,.233,'4 癌平均上调，3 癌显著；约 73%–79% 患者同向。',13,G,True)
text(.447,.201,'* 名义 P < 0.05；PDAC P = 0.0777，未达显著。',10.2,M)
line(.05,.95,.176)
text(.05,.131,'由此提出问题',13,G,True)
text(.183,.131,'这种跨癌重复，是否对应恶性上皮细胞中的共同上调？',19,G,True)
text(.05,.071,'这是进一步关注的候选，不代表最优靶点；ccRCC 已见反向，因此不称“泛癌一致上调”。',10.5,M)
text(.05,.039,'右侧为同一 CAMP 队列内的表达支持，不是独立验证或机制证明；下一页展示四癌患者配对细节。',9.8,M)
text(.95,.041,'09',14,G,True,ha='right')
f.canvas.draw();rend=f.canvas.get_renderer()
for t in f.texts:
 bb=t.get_window_extent(rend);assert bb.x0>=0 and bb.y0>=0 and bb.x1<=f.bbox.width and bb.y1<=f.bbox.height,t.get_text()
name='09_为什么关注LYPLA1'
for ext in ['png','pdf','svg']:
 p=OUT/f'{name}.{ext}';f.savefig(p,dpi=200)
 if ext=='svg':p.write_text('\n'.join(s.rstrip() for s in p.read_text(encoding='utf8').splitlines())+'\n',encoding='utf8')
 shutil.copy2(p,VIEW/p.name)
plt.close(f)
spec=dict(scope='LYPLA1 motivation page after six-cancer overview; four cancers shown',logic=['CAMP metabolite discovery','direct biochemical candidate mapping','four-cancer expression and patient direction motivates LYPLA1 follow-up','question of malignant epithelial source'],selection_claim='Not claimed as unique or best candidate; four-cancer scope does not erase ccRCC reverse result',statistics='Frozen paired RNA results reused; n and P crosschecked against page8; no new tests or patient filtering',nominal_p=.05,paired_detail='Move existing paired detail to next page10',software={'pandas':pd.__version__,'matplotlib':matplotlib.__version__})
(OUT/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf8')
(OUT/'README_CN.md').write_text('''# 为什么关注LYPLA1：专题引入页

本轮问题：用户要求先讲清为什么关注到LYPLA1，再展开验证结果。本页作为第9页，既有四癌患者配对连线图作为第10页。

输入与范围：冻结的四癌配对RNA汇总与第8页患者方向计数。只读取既有公开汇总，原效应、P、入选患者数不变；无新推断检验。

实际结果：BRCA34/45（75.6%）、COAD26/33（78.8%）、PDAC8/11（72.7%）、PRAD33/43（76.7%）位患者上调。四癌平均表达变化为正，三癌名义P<0.05；PDAC P=0.0777。右侧P与n逐项核对第8页。

新手解释：左侧解释候选来自代谢异常及直接关系映射，右侧展示表达和患者同向线索，底部提出恶性上皮细胞来源问题。这是选择进一步研究对象的理由，不是将表达相关升级为机制验证。

限制/反证：LYPLA1不是唯一并列候选，没有声称其统计排名优于其他基因。ccRCC反向仍在概览页与本页脚注提示，不能因四癌展示就称泛癌一致。配对患者上调比例不等于每位患者各自显著。四癌未合并患者数据或P值；同一CAMP队列内的支持不是独立验证。

当前决定：生成图片、PDF和SVG，不生成PPT；承接第8页，后接配对细节，再进入单细胞。

复现：python code/ppt_page9_focus_reason_v1.py。
''',encoding='utf8')
(OUT/'validation.json').write_text(json.dumps(dict(status='DONE',original_n_and_p_crosschecked=True,four_positive_means=True,three_nominal_significant=True,patient_percent_range=[float(d.up_fraction.min()),float(d.up_fraction.max())],new_tests=0,visual_review='PENDING'),indent=2),encoding='utf8')
for n in ['README_CN.md','four_cancer_evidence.tsv']:shutil.copy2(OUT/n,VIEW/n)
(VIEW/'index.html').write_text(f'<!doctype html><meta charset="utf-8"><title>CAMP · 为什么关注LYPLA1</title><style>body{{background:#e8eeea;margin:0;font-family:system-ui}}main{{max-width:1440px;margin:24px auto}}img{{width:100%;display:block;box-shadow:0 8px 25px #163e3320}}a{{color:#00553b}}</style><main><h2>第9页：为什么关注 LYPLA1？</h2><img src="{name}.png"><p><a href="{name}.pdf">矢量 PDF</a> · <a href="four_cancer_evidence.tsv">四癌数据</a></p></main>',encoding='utf8')
print(VIEW)
