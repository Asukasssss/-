"""Concise evidence-convergence slide; all numerical evidence remains in v1."""
from pathlib import Path
import json,shutil,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle,FancyArrowPatch
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/BRCA/07_INTEGRATION'
RUN='20260929T_lypla1_evidence_summary_v2';OUT=BASE/RUN;OUT.mkdir(exist_ok=True)
OLD=BASE/'20260929T_lypla1_evidence_summary_v1'
VIEW=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_LYPLA1_evidence_summary_v2');VIEW.mkdir(exist_ok=True)
for n in ['msyh.ttc','msyhbd.ttc']:font_manager.fontManager.addfont('C:/Windows/Fonts/'+n)
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42})
G='#00553B';INK='#183D34';M='#6A7973';LINE='#DDE7E0'
f=plt.figure(figsize=(16,9),facecolor='white')
def t(x,y,s,size=16,col=INK,bold=False,ha='left'):return f.text(x,y,s,fontsize=size,color=col,weight='bold' if bold else 'normal',va='center',ha=ha,linespacing=1.7)
def line(x1,x2,y,col=LINE,w=.8):f.add_artist(Line2D([x1,x2],[y,y],transform=f.transFigure,color=col,lw=w))
ax=f.add_axes([.044,.882,.168,.10]);ax.imshow(plt.imread(BASE/'20260929T170000Z_page8_genes_six_v2/sysu_logo.png'));ax.axis('off')
t(.95,.93,'CAMP  /  组会汇报',10,M,ha='right');line(.05,.95,.883);line(.05,.103,.883,G,2)
t(.05,.821,'LYPLA1：从代谢候选到多层表达证据',27,G,True)
rows=[('01','代谢发现','四癌异常代谢物经关系映射，均指向 LYPLA1'),('02','组织表达','四癌平均上调，其中 BRCA、COAD、PRAD 显著'),('03','患者重复','四癌上皮比较中，多数患者呈同向升高'),('04','空间定位','BRCA、PRAD 代表切片提供癌区定位补充')]
ys=[.684,.541,.398,.255]
for (num,title,body),y in zip(rows,ys):
 t(.06,y,num,23,'#A5BDAF',True);t(.107,y,title,20,G,True)
 t(.107,y-.054,body,14.2,INK)
 if num!='04':line(.107,.65,y-.103)
# Lines connect evidence summaries, not biochemical pathways.
f.add_artist(Line2D([.687,.707,.707,.687],[.699,.699,.230,.230],transform=f.transFigure,color='#91AD9C',lw=1.5))
f.add_artist(FancyArrowPatch((.707,.465),(.744,.465),transform=f.transFigure,arrowstyle='-|>',mutation_scale=18,color='#91AD9C',lw=1.5))
f.add_artist(Rectangle((.755,.250),.195,.430,transform=f.transFigure,facecolor=G,edgecolor='none',zorder=-1))
t(.8525,.597,'LYPLA1',30,'white',True,ha='center')
line(.785,.920,.543,'#679A85',1)
t(.8525,.468,'跨癌反复出现\n具有上皮表达支持',16,'white',ha='center')
t(.8525,.334,'值得继续关注的候选',14,'#E0EEE5',True,ha='center')
line(.05,.95,.133)
t(.05,.086,'现阶段：多层证据支持候选价值，共同代谢机制尚未建立。',16,G,True)
t(.95,.038,'证据汇总',10,G,True,ha='right')
stem='LYPLA1_多层证据收拢'
for ext in ['png','pdf','svg']:f.savefig(OUT/(stem+'.'+ext),dpi=180,facecolor='white');shutil.copy2(OUT/(stem+'.'+ext),VIEW/(stem+'.'+ext))
plt.close(f)
for name in ['source_metabolite.tsv','source_rna.tsv','source_sc.tsv','source_spatial.json']:shutil.copy2(OLD/name,OUT/name)
shutil.copy2(OLD/'README_CN.md',OUT/'evidence_details_CN.md')
(OUT/'source_manifest.tsv').write_text('source\tsha256\n'+''.join(f'{p.relative_to(ROOT)}\t{hashlib.sha256(p.read_bytes()).hexdigest()}\n' for p in sorted(OLD.glob('source_*'))),encoding='utf-8')
(OUT/'analysis_spec.json').write_text(json.dumps({'version':'v2','change':'Replace detailed table with four evidence statements converging on candidate conclusion','statistics':'unchanged; none recalculated','scope':['BRCA','COAD','PDAC','PRAD'],'software':matplotlib.__version__},ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'README_CN.md').write_text('''# LYPLA1 证据收拢页 v2
## 本轮问题
用户要求减少小字，收拢证据，不再重复逐癌细节表。
## 输入与范围
复用 v1 汇总表；只改展示，不改分析。四癌范围保持 BRCA、COAD、PDAC、PRAD。
## 实际结果
四条大字证据向 LYPLA1 候选结论汇合，输出 PNG/PDF/SVG。
## 新手解释
代谢关系映射发现候选；组织 RNA 显示同向趋势；同患者上皮比较提供重复支持；空间图补充定位。
## 限制/反证
PDAC 组织 RNA 未显著，因此仅写其中三癌显著。上皮参照、描述性检验范围、PDAC 阳性核敏感性及空间患者限制完整保存在 evidence_details_CN.md；未将这些证据改称共同机制或正式泛癌验证。
## 当前决定
展示页精简，详细统计与限制保留在前页和附属文档。未生成 PPTX。
## 下一步
组会用此页收拢 LYPLA1 证据。
## 复现命令
python code/ppt_lypla1_evidence_summary_v2.py
''',encoding='utf-8')
(OUT/'validation.json').write_text(json.dumps({'status':'PARTIAL','render':'PASS','visual_review':'pending','new_tests':False},indent=2),encoding='utf-8')
html=f'''<!doctype html><meta charset="utf-8"><title>LYPLA1 证据收拢</title><style>body{{margin:0;background:#e9eeeb;color:#183d34;font-family:"Microsoft YaHei",sans-serif}}nav{{padding:18px 4%;background:white}}main{{max-width:1680px;margin:24px auto}}img{{width:100%;display:block}}a{{color:#00553b}}</style><nav>LYPLA1 · 证据收拢　<a href="{stem}.pdf">PDF</a>　<a href="{stem}.png">高清 PNG</a>　<a href="{stem}.svg">SVG</a>　<a href="evidence_details_CN.md">详细证据</a></nav><main><img src="{stem}.png"></main>'''
for d in [OUT,VIEW]:(d/'index.html').write_text(html,encoding='utf-8')
shutil.copy2(OUT/'evidence_details_CN.md',VIEW/'evidence_details_CN.md')
print(VIEW)
