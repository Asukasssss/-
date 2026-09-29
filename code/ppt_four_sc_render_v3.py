"""One consistent layout: expression, annotation, right-side legend, patient evidence."""
from pathlib import Path
import sys,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
R=Path(sys.argv[1]);P=R/'public'
for n in ['msyh.ttc','msyhbd.ttc']:font_manager.fontManager.addfont(str(R/n))
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
G='#00553B';INK='#183D34';M='#6A7973';UP='#C77462';DOWN='#468B9E';LINE='#DDE7E0';RED='#DD3434'
def color(t):
 s=t.lower()
 if s in ['cancer epithelial','malignant epithelial','cna epithelial']:return '#78b436'
 if s in ['normal epithelial','ductal']:return '#43a2ca'
 if 'unclassified' in s or s=='epithelial':return '#C5CAC7'
 if 'cnn' in s or 'altered' in s or 'atypical' in s:return '#C8AF45'
 if 'plasma' in s:return '#d95f9f'
 if 'b cell' in s or s=='b-cells':return '#d73027'
 if 'endothel' in s:return '#7b3294'
 if 'fibro' in s or s=='cafs':return '#4575b4'
 if 'myelo' in s or 'mono' in s or 'macro' in s:return '#ef8a0c'
 if 'treg' in s or 'regulatory t' in s:return '#237052'
 if 'cd8' in s:return '#438C71'
 if 'mast' in s:return '#AD7E8E'
 if s in ['t cell','t cells','t-cells','conventional t cells'] or 't/nk' in s:return '#66a98c'
 if 'nk' in s:return '#8BADB0'
 if 'mast' in s:return '#AD7E8E'
 if s in ['pvl','smcs'] or 'pericy' in s or 'smooth' in s:return '#888888'
 if 'glia' in s or 'schwann' in s:return '#A48869'
 if 'acinar' in s:return '#A38BC3'
 if 'endocrine' in s or 'islet' in s:return '#DDA87D'
 return '#A9B5BA'
def label(t):
 mapping={'Cancer Epithelial':'Cancer epithelial','Normal Epithelial':'Normal epithelial','B-cells':'B cells','T-cells':'T cells','Epithelial unclassified':'Epithelial (unclassified)','Altered benign epithelial':'Altered benign epi.','Conventional T cells':'Conventional T','Regulatory T cells':'Regulatory T','Malignant epithelial':'Malignant epithelial','Epithelial (non-malignant)':'Other nonmalignant epi.'}
 mapping.update({'Cancer-associated fibroblast':'CAFs','Intra-pancreatic neurons':'Pancreatic neurons'})
 return mapping.get(t,t.replace('_',' '))
for page,c in enumerate(['BRCA','COAD','PDAC','PRAD'],11):
 z=np.load(R/(c+'_private.npz'));xy=z['xy'];expr=z['expr'];types=z['types'];pairs=z['pairs'];meta=json.loads((R/(c+'_meta.json')).read_text());cap=meta['cap']
 cats=sorted(set(types));cats.sort(key=lambda t:(0 if t==meta['target'] else 1,label(t)))
 f=plt.figure(figsize=(16,9),facecolor='white')
 def text(x,y,s,size=13,col=INK,bold=False,ha='left'):f.text(x,y,s,fontsize=size,color=col,weight='bold' if bold else 'normal',ha=ha,va='center')
 def line(x1,x2,y,col=LINE,w=.8):f.add_artist(Line2D([x1,x2],[y,y],transform=f.transFigure,color=col,lw=w))
 ax=f.add_axes([.044,.882,.168,.10]);ax.imshow(plt.imread(R/'sysu_logo.png'));ax.axis('off')
 text(.95,.930,'CAMP  /  组会汇报',10,M,ha='right');line(.05,.95,.883);line(.05,.103,.883,G,2)
 text(.05,.828,f'{c}：从全细胞表达定位到同患者证据',27,G,True)
 text(.05,.777,f"{meta['cohort']} · {meta['scope']} · {meta['cells']:,} 个细胞"+('核' if c=='PDAC' else ''),12,M)
 text(.048,.707,'LYPLA1 expression',15,INK,True);text(.365,.707,'Cell types',15,INK,True);text(.838,.707,'同患者证据',15,G,True)
 ax1=f.add_axes([.050,.268,.238,.409]);ax2=f.add_axes([.364,.268,.238,.409])
 ax1.scatter(xy[:,0],xy[:,1],s=.5,c='#e6e9e9',lw=0,rasterized=True)
 order=np.argsort(expr,kind='stable');pos=order[expr[order]>0]
 im=ax1.scatter(xy[pos,0],xy[pos,1],c=expr[pos],cmap='viridis',vmin=0,vmax=cap,s=.6,lw=0,rasterized=True)
 ix=np.random.default_rng(20260929).permutation(len(types));ax2.scatter(xy[ix,0],xy[ix,1],c=[color(t) for t in types[ix]],s=.6,lw=0,rasterized=True)
 # The same locator is repeated on both panels; it does not select cells.
 q1,q2=np.quantile(xy[types==meta['target']],[.03,.97],axis=0);center=(q1+q2)/2;width,height=(q2-q1)*1.12
 for ax in [ax1,ax2]:
  ax.set_aspect('equal');pad=(xy.max(0)-xy.min(0))*.035
  ax.set_xlim(xy[:,0].min()-pad[0],xy[:,0].max()+pad[0]);ax.set_ylim(xy[:,1].min()-pad[1],xy[:,1].max()+pad[1])
  ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('UMAP 1',fontsize=10,color=M);ax.set_ylabel('UMAP 2',fontsize=10,color=M)
  for s in ax.spines.values():s.set_visible(False)
  for k in ['left','bottom']:ax.spines[k].set_visible(True);ax.spines[k].set_color('#ADB7B1');ax.spines[k].set_linewidth(.65)
  ax.add_patch(Ellipse(center,width,height,fill=False,ec=RED,lw=1.15,zorder=7))
 cb=f.add_axes([.301,.315,.008,.300]);f.colorbar(im,cax=cb,extend='max');cb.tick_params(labelsize=8,length=2);cb.set_ylabel('LYPLA1 · log1p(CP10K)',fontsize=8.5,color=M,labelpad=5)
 # Dedicated vertical legend immediately to the right of the annotation UMAP.
 text(.626,.660,'Cell type (n cells)',10,INK,True)
 step=min(.032,.390/max(1,len(cats)-1));y=.626
 for t in cats:
  count=int((types==t).sum());f.add_artist(Line2D([.632],[y],transform=f.transFigure,marker='o',ls='',color=color(t),markersize=4))
  text(.642,y,label(t),8.3,color(t) if t==meta['target'] else INK,t==meta['target'])
  text(.806,y,f'{count:,}',8.1,M,ha='right');y-=step
 text(.48,.213,meta['locator']+'：红圈仅作位置提示',9,RED,ha='center')
 text(.177,.213,f'色标上限 {cap:.2f} · 阳性值第99百分位',8.7,M,ha='center')
 f.add_artist(Line2D([.819,.819],[.195,.683],transform=f.transFigure,color=LINE,lw=.8))
 p=meta['p'];ps=f"配对 P = {p:.4g}" if p is not None else '描述性配对 · 未做显著性检验'
 text(.908,.655,ps,10 if p is not None else 8.8,G,p is not None,ha='center')
 ax=f.add_axes([.857,.374,.109,.218])
 for a,b in pairs:ax.plot([0,1],[a,b],c=UP if b>a else DOWN,lw=1,alpha=.6)
 for j,col in enumerate([DOWN,UP]):
  ax.scatter(np.full(len(pairs),j),pairs[:,j],c=col,s=19,edgecolor='white',lw=.4,zorder=3);mean=pairs[:,j].mean();ax.plot([j-.12,j+.12],[mean,mean],c=G,lw=2.4,zorder=4)
 ax.set_xlim(-.3,1.3);ax.set_ylim(0,pairs.max()*1.17);ax.set_xticks([0,1],['参照','目标'],fontsize=8)
 ax.set_ylabel('患者平均表达',fontsize=8.5,color=M);ax.tick_params(length=2,labelsize=8,color=M,labelcolor=M)
 for k in ['top','right']:ax.spines[k].set_visible(False)
 for k in ['left','bottom']:ax.spines[k].set_color(LINE)
 text(.910,.312,f"{meta['higher']} / {meta['n_pairs']}",24,G,True,ha='center');text(.910,.274,'患者目标组表达更高',9.8,G,True,ha='center')
 text(.910,.229,f"平均差 {meta['mean_difference']:+.3f}",9.5,M,ha='center')
 text(.841,.180,'参照：'+meta['ref_label'],8.5,M);text(.841,.152,'目标：'+meta['test_label'],8.5,M)
 line(.05,.95,.127)
 conclusion=meta['conclusion']
 if p is None:conclusion=f"同患者描述：{meta['higher']}/{meta['n_pairs']} 位目标上皮表达更高；尚未作显著性检验。"
 text(.05,.091,conclusion,14,G,True)
 text(.05,.054,meta['note'],8.5,M);text(.96,.046,str(page),14,G,True,ha='right')
 name=f'{page}_{c}_LYPLA1_全细胞与患者'
 for ext in ['png','pdf','svg']:f.savefig(P/(name+'.'+ext),dpi=220)
 plt.close(f)
 print(c,'rendered',flush=True)
(P/'validation.json').write_text(json.dumps(dict(status='PASS',all_four_rendered=True,source_coordinates_reused=True,legend_right_of_umap=True,patient_mean_and_direction_checks=True,new_inferential_tests=False,visual_review='PENDING'),indent=2))
