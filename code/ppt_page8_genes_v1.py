"""Read-only ranking of current candidate RNA results across six cancers."""
from pathlib import Path
import hashlib, io, json, shutil, subprocess
import pandas as pd
import numpy as np
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--scope',choices=['four','six'],default='four');args=parser.parse_args()
ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'results/BRCA/07_INTEGRATION/20260928T150000Z_ppt_redraw_v1'
OUT=ROOT/'results/BRCA/07_INTEGRATION/20260929T160000Z_page8_genes_v1'
VIEW=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_page08_genes_v1')
if args.scope=='six':
 OUT=ROOT/'results/BRCA/07_INTEGRATION/20260929T170000Z_page8_genes_six_v2'
 VIEW=VIEW.parent/'CAMP_page08_genes_six_v2'
OUT.mkdir(parents=True,exist_ok=True);VIEW.mkdir(parents=True,exist_ok=True)
manifest=pd.read_csv(OLD/'source_manifest.tsv',sep='\t')
tables={};sources=[]
for c in ['brca','coad','pdac','prad']:
 r=manifest[manifest.id==c+'_gene'].iloc[0].to_dict()
 b=(OLD/'sources'/f'{c}_gene.tsv').read_bytes()
 assert hashlib.sha256(b).hexdigest()==r['sha256']
 tables[c]=pd.read_csv(io.BytesIO(b),sep='\t');sources.append(r)
for c,commit,path in [('gbm','2174aa965e29cc0fbb10925bc2c217198960829e','results/GBM/07_INTEGRATION/20260925T145000Z_integration_v1/candidate_genes_integrated.tsv'),('ccrcc','02593f745fc966c975ab1a25e3e7ce74f770efc6','results/ccRCC/07_INTEGRATION/20260925T163000Z_integration_v1/candidate_genes_integrated.tsv')]:
 b=subprocess.check_output(['git','show',f'{commit}:{path}'],cwd=ROOT)
 tables[c]=pd.read_csv(io.BytesIO(b),sep='\t')
 sources.append(dict(id=c+'_gene',git_commit=commit,git_path=path,sha256=hashlib.sha256(b).hexdigest(),bytes=len(b),url=f'https://github.com/Asukasssss/-/blob/{commit}/{path}'))
pd.DataFrame(sources).to_csv(OUT/'source_manifest.tsv',sep='\t',index=False)
config=[('BRCA','brca','v2_in_current190_direct_pool','v2_paired_RNA_effect','v2_paired_RNA_p_value','v2_paired_RNA_status','v2_paired_RNA_n'),('COAD','coad','current_direct','RNA_effect','RNA_p_value','RNA_status','RNA_n'),('GBM','gbm','current_pool','RNA_effect','RNA_p_value','RNA_status','RNA_n'),('PDAC','pdac','in_current_pool','RNA_effect','RNA_p','RNA_status','RNA_n'),('PRAD','prad','current_pool','RNA_effect','RNA_p','RNA_status','RNA_n')]
for c in ['ccRCC3','ccRCC4']:config.append((c,'ccrcc','current_in_'+c,c+'_RNA_effect',c+'_RNA_p_value',c+'_RNA_status',c+'_RNA_n'))
records=[];counts=[]
for cohort,tab,flag,eff,pval,status,n in config:
 d=tables[tab];d=d[d[flag].eq(True)].copy();assert not d.gene.duplicated().any()
 counts.append(dict(cohort=cohort,current_pool=len(d),done=int(d[status].eq('DONE').sum())))
 for _,r in d.iterrows():
  good=r[status]=='DONE' and np.isfinite(r[eff]) and np.isfinite(r[pval])
  records.append(dict(gene=r.gene,cohort=cohort,cancer='ccRCC' if cohort.startswith('ccRCC') else cohort,effect=r[eff],p=r[pval],n=r[n],source_status=r[status],status='EVALUABLE' if good else 'NOT_EVALUABLE'))
d=pd.DataFrame(records);assert d.groupby(['gene','cohort']).size().max()==1
order=[x[0] for x in config]
d=d.set_index(['gene','cohort']).reindex(pd.MultiIndex.from_product([sorted(d.gene.unique()),order],names=['gene','cohort'])).reset_index()
d['status']=d.status.fillna('OUTSIDE_CURRENT_POOL');d['cancer']=d.cohort.replace({'ccRCC3':'ccRCC','ccRCC4':'ccRCC'})
if args.scope=='four':
 order=['BRCA','COAD','PDAC','PRAD'];d=d[d.cohort.isin(order)].copy()
rows=[]
for g,z in d.groupby('gene'):
 z=z[z.status=='EVALUABLE'];up=z[(z.p<.05)&(z.effect>0)];dn=z[(z.p<.05)&(z.effect<0)];u=up.cancer.nunique();v=dn.cancer.nunique();direction='up' if u>v else 'down' if v>u else 'tie'
 rows.append(dict(gene=g,same_significant_cancers=max(u,v),opposite_significant_cancers=min(u,v),dominant_direction=direction,evaluable_cancers=z.cancer.nunique(),same_significant_cohorts=len(up if direction=='up' else dn),evaluable_cohorts=len(z)))
ranking=pd.DataFrame(rows).sort_values(['same_significant_cancers','opposite_significant_cancers','evaluable_cancers','same_significant_cohorts','gene'],ascending=[False,True,False,False,True])
sel=ranking[(ranking.dominant_direction!='tie')&(ranking.same_significant_cancers==ranking.same_significant_cancers.max())].copy()
display=d[d.gene.isin(sel.gene)].copy()
for name,a in [('ranking_all',ranking),('display_genes',sel),('display_cells',display),('all_cells',d),('coverage',pd.DataFrame(counts))]:a.to_csv(OUT/f'{name}.tsv',sep='\t',index=False)
print(pd.DataFrame(counts).to_string(index=False));print(ranking.head(25).to_string(index=False));print(display.pivot(index='gene',columns='cohort',values='effect').to_string())

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle
for fp in ['C:/Windows/Fonts/msyh.ttc','C:/Windows/Fonts/msyhbd.ttc']:font_manager.fontManager.addfont(fp)
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
G='#00553B';INK='#183D34';M='#6A7973';UP='#C77462';DOWN='#468B9E';LINE='#DDE7E0'
shutil.copy2(ROOT/'results/BRCA/07_INTEGRATION/20260929T150000Z_page7_same_direction_v3/sysu_logo.png',OUT/'sysu_logo.png')
f=plt.figure(figsize=(16,9),facecolor='white')
def text(x,y,s,size=15,col=INK,bold=False,ha='left'):return f.text(x,y,s,fontsize=size,color=col,weight='bold' if bold else 'normal',ha=ha,va='center')
def line(x1,x2,y,c=LINE,w=.8):f.add_artist(plt.Line2D([x1,x2],[y,y],transform=f.transFigure,color=c,lw=w))
ax=f.add_axes([.044,.874,.177,.106]);ax.imshow(plt.imread(OUT/'sysu_logo.png'));ax.axis('off')
text(.95,.928,'CAMP  /  组会汇报',10,M,ha='right');line(.05,.95,.876);line(.05,.103,.876,G,2.2)
text(.05,.813,'哪些候选基因在更多癌种中同向改变？',29,G,True)
text(.052,.755,'RNA 肿瘤—参照比较  ·  展示同向 P < 0.05 癌种数最多的全部并列基因',13.5,M)
ax=f.add_axes([.05,.23,.90,.47]);ax.set_xlim(-2.5,9.05);ax.set_ylim(len(sel)-.15,-1.2);ax.axis('off')
xx=np.linspace(.1,6.1,len(order))
for j,c in enumerate(order):ax.text(xx[j],-.65,c,ha='center',va='center',fontsize=12,color=G,weight='bold')
ax.text(-2.4,-.65,'候选基因',fontsize=12,color=G,va='center',weight='bold')
ax.text(8.05,-.76,'同向显著 / 可评估',fontsize=10.7,ha='center',va='center',color=G,weight='bold');ax.text(8.05,-.37,'癌种数',fontsize=10,ha='center',va='center',color=M)
ax.plot([-2.5,9.05],[-.12,-.12],color=LINE,lw=.8)
for j in range(len(order)):ax.plot([xx[j],xx[j]],[0,len(sel)-.2],c='#ECF0EC',lw=.6,zorder=-3)
for i,r in enumerate(sel.itertuples()):
 y=i+.3
 if i%2==0:ax.add_patch(Rectangle((-2.5,y-.43),11.55,.86,fc='#F7F9F7',ec='none',zorder=-4))
 ax.text(-2.4,y,r.gene,fontsize=15,color=INK,va='center',weight='bold')
 for j,c in enumerate(order):
  a=display[(display.gene==r.gene)&(display.cohort==c)].iloc[0]
  if a.status=='EVALUABLE':
   color=UP if a.effect>0 else DOWN;ax.scatter(xx[j],y,s=175,facecolors=color if a.p<.05 else 'white',edgecolors=color,linewidths=1.5,zorder=3)
  else:ax.text(xx[j],y,'—' if a.status=='OUTSIDE_CURRENT_POOL' else '×',ha='center',va='center',color='#A5B5AB',fontsize=18)
 color=UP if r.dominant_direction=='up' else DOWN
 ax.text(7.35,y,'↑' if r.dominant_direction=='up' else '↓',fontsize=19,color=color,ha='center',va='center')
 ax.text(8.05,y,f'{r.same_significant_cancers} / {r.evaluable_cancers}',fontsize=20,color=color,ha='center',va='center',weight='bold')
ax.plot([6.9,6.9],[-1,len(sel)-.2],color=LINE,lw=.8)
legend=f.add_axes([.14,.175,.80,.035]);legend.set_xlim(0,12);legend.set_ylim(0,1);legend.axis('off')
for x,c,l in [(0,UP,'升高'),(1.5,DOWN,'降低')]:legend.scatter(x+.1,.5,s=85,c=c);legend.text(x+.32,.5,l,fontsize=10.5,color=M,va='center')
legend.scatter(3.05,.5,s=85,c=M);legend.text(3.28,.5,'实心 P < 0.05',fontsize=10.5,color=M,va='center')
legend.scatter(5.7,.5,s=85,facecolors='white',edgecolors=M,linewidths=1.3);legend.text(5.93,.5,'空心 P ≥ 0.05',fontsize=10.5,color=M,va='center')
legend.text(8.55,.5,'—',fontsize=15,color='#A5B5AB',va='center');legend.text(8.9,.5,'当前候选池未覆盖',fontsize=10.5,color=M,va='center')
line(.05,.95,.15)
if args.scope=='six':
 text(.05,.113,'同向显著最多：4 个癌种',14,G,True)
 text(.40,.113,'6 个基因并列；不同癌种仍可能出现反向结果。',12,INK)
 text(.05,.077,'LYPLA1 在 ccRCC3、SLC7A11 在 GBM 显著降低；同向计数不代表所有癌种一致。',10.5,M)
 text(.05,.041,'范围：当前直接关系候选池，非全转录组。ccRCC 两队列计一个癌种；GBM 为非配对跨来源参照。',9.5,M)
else:
 text(.05,.113,'同向显著最多：3 个癌种',14,G,True)
 text(.40,.113,'7 个基因并列；没有候选达到四癌同向显著。',12,INK)
 text(.05,.077,'LYPLA1 四癌均升高、三癌显著；SLC29A2 在 PDAC 方向相反但未显著。',10.5,M)
 text(.05,.041,'范围：四癌当前直接关系候选池，非全转录组。AHCY、SLC6A6 各有一个显著反向癌种。',9.5,M)
text(.95,.041,'08',14,G,True,ha='right')
f.canvas.draw();rend=f.canvas.get_renderer()
for t in f.texts:
 b=t.get_window_extent(rend);assert b.x0>=0 and b.y0>=0 and b.x1<=f.bbox.width and b.y1<=f.bbox.height
name='08_候选基因同向最多'
for ext in ['png','pdf','svg']:
 p=OUT/f'{name}.{ext}';f.savefig(p,dpi=200)
 if ext=='svg':p.write_text('\n'.join(x.rstrip() for x in p.read_text(encoding='utf8').splitlines())+'\n',encoding='utf8')
 shutil.copy2(p,VIEW/p.name)
plt.close(f)
spec=dict(scope=args.scope,cohorts=order,selection='all genes tied at maximum same-direction nominal P<0.05 cancer count',ties='fewer opposite significant cancers; more evaluable cancers; more same significant cohorts; symbol alphabetically',pool='current direct biochemical candidates per cancer only; not whole transcriptome',p=.05,new_patient_tests=False,ccRCC_count='one cancer per direction; either cohort can contribute; two significant opposing cohorts contribute to both directional counts, never add as distinct cancers',gene_identity='exact curated gene symbols in integration tables; assert unique per cohort; no new synonym remapping',software=dict(pandas=pd.__version__,numpy=np.__version__,matplotlib=matplotlib.__version__),source_fields=config)
(OUT/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf8')
assert len(display)==len(sel)*len(order)
assert ((d.loc[d.status=='EVALUABLE','p']>=0)&(d.loc[d.status=='EVALUABLE','p']<=1)).all()
(OUT/'validation.json').write_text(json.dumps(dict(status='DONE',unique_gene_cohort=True,source_hashes_verified=True,all_maximum_ties_displayed=True,selected_genes=sel.gene.tolist(),selected_count=len(sel),display_cells=len(display),scope=args.scope,visual_review='See rendered image; no patient-level retesting'),indent=2),encoding='utf8')
(OUT/'README_CN.md').write_text(f'''# 第8页：候选基因同向重复

本轮问题：找出RNA肿瘤—参照比较中，同向名义P<0.05癌种数最多的候选基因。

输入与范围：{args.scope}；各癌种当前直接关系候选池，不是全转录组。历史、条件性、未解决关系不入池。源表版本及SHA256见source_manifest.tsv；字段映射见analysis_spec.json。原效应、P值及样本数直接复用，不重算患者统计。

实际结果：最高同向显著癌种数为{int(sel.same_significant_cancers.max())}；并列基因为{', '.join(sel.gene)}。全部并列展示。并列内部依次按反向显著癌种更少、可评估癌种更多、同向显著队列更多、符号字母排序；此内部顺序不是效应强度排名。

新手解释：暖色升高、蓝色降低；实心P<0.05，空心P>=0.05。横线表示该癌种当前候选池未覆盖，不代表基因未表达、无差异或缺RNA数据。×表示候选在池中但本轮不能评估。分母是当前池内可评估癌种数。ccRCC两队列分别画，但每方向只计一个癌种；若两队列各有显著反向结果，方向计数分别记录、不能相加当不同癌种。统计RNA差异，不混入代谢物RNA相关或单细胞结果。

限制/反证：不同癌种候选池和发现覆盖不同，不能将当前池外当阴性。效应尺度沿用各源分析，点大小不编码效应，不跨癌混算原效应。GBM为肿瘤与GTEx非配对比较，存在来源混杂。未使用FDR作为展示门槛；名义显著仅作探索性证据。反向与不显著结果均保留。

当前决定：制作图片及矢量文件，不生成PPT；后续LYPLA1专题仍按此前约定只比较BRCA、COAD、PDAC、PRAD。

下一步：与代谢物同向结果衔接，再进入LYPLA1专题。

复现：python code/ppt_page8_genes_v1.py --scope {args.scope}
''',encoding='utf8')
for n in ['display_genes.tsv','display_cells.tsv','ranking_all.tsv','analysis_spec.json','README_CN.md']:shutil.copy2(OUT/n,VIEW/n)
(VIEW/'index.html').write_text(f'<!doctype html><meta charset="utf-8"><title>CAMP 第8页 · 候选基因</title><style>body{{background:#e8eeea;margin:0;font-family:system-ui}}main{{max-width:1440px;margin:24px auto}}img{{width:100%;display:block;box-shadow:0 8px 25px #163e3320}}a{{color:#00553b}}</style><main><h2>第8页：同向显著癌种数最多的候选基因</h2><img src="{name}.png"><p><a href="{name}.pdf">矢量 PDF</a> · <a href="display_genes.tsv">排名与覆盖</a> · <a href="README_CN.md">统计范围</a></p></main>',encoding='utf8')
files=[p for p in OUT.iterdir() if p.is_file() and p.name!='checksums.tsv']
pd.DataFrame([dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]).to_csv(OUT/'checksums.tsv',sep='\t',index=False)
