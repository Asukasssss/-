"""Server-side redraw of four user-selected sections; no statistical recomputation."""
from pathlib import Path
import sys,json,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from PIL import Image
R=Path(sys.argv[1]);P=R/'public';P.mkdir(exist_ok=True)
B=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716');C=B/'results/collaborative'
AS=C/'BRCA/A/20260929T071500Z_four_sc_v3'
for n in ['msyh.ttc','msyhbd.ttc']:font_manager.fontManager.addfont(str(AS/n))
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':False})
G='#00553B';INK='#183D34';M='#6A7973';LINE='#DDE7E0';ZERO='#E5E8EC';CAP=3.0
COL={'Tumor':'#C94B54','Cancer':'#C94B54','Benign':'#3276A5','Immune':'#D1AE35','DCIS':'#836BB1','Unannotated':'#BFC6C9','Mixed':'#E4E8E8','Stroma':'#BFC6C9','PIN':'#A98C24','Transition_State':'#258F88','Vessel':'#258F88','Necrosis':'#777777'}
LABEL={'Tumor':'癌区','Cancer':'癌区','Benign':'良性腺体','Immune':'免疫区','DCIS':'原位癌','Unannotated':'未标注区','Mixed':'混合区','Stroma':'间质','PIN':'PIN','Transition_State':'过渡区','Vessel':'血管','Necrosis':'坏死'}
COL['Fat']='#C4AE94';LABEL['Fat']='脂肪区'
COL['Nerve']='#749487';LABEL['Nerve']='神经区'
sources=[];checks=[];data=[]
def source(p):
 sources.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
for s in ['V19T26-012_B1','V19T26-032_B1']:
 r=C/'BRCA/A/20260928T060000Z_cta_lypla1_v1';p=r/(s+'_spots_SERVER_ONLY.tsv.gz');source(p);d=pd.read_csv(p,sep='\t');d=d[d.qc.astype(str).str.lower().isin(['true','1'])]
 ip=B/'data/candidates/BRCA_CTA2025/spaceranger_output'/s/'outs/spatial/tissue_lowres_image.png';source(ip)
 data.append(dict(cancer='BRCA',sample=s,patient=d.patient.iloc[0],img=np.asarray(Image.open(ip)),xy=d[['x','y']].to_numpy(),expr=d.lognorm.to_numpy(),region=d.region.to_numpy()))
r=C/'PRAD/B/20260928T101700Z_lypla1_erickson_v1';p=r/'private/spot_values.tsv.gz';source(p);d=pd.read_csv(p,sep='\t')
for s in ['H2_5','H2_1']:
 x=d[d.section==s];folder=r/'source/Patient1'/s;sp=folder/'scalefactors_json.json';source(sp);sf=json.loads(sp.read_text());ip=folder/(s+'_tissue_hires_image.png');source(ip)
 data.append(dict(cancer='PRAD',sample=s,patient='Patient 1',img=np.asarray(Image.open(ip)),xy=x[['pixel_col','pixel_row']].to_numpy()*sf['tissue_hires_scalef'],expr=x.expression.to_numpy(),region=x.region.to_numpy()))
for e in data:
 for group in sorted(set(e['region'])):
  if group not in COL:COL[group]='#A3B7AB';LABEL[group]=str(group)
print('Observed categories:',sorted(set(g for e in data for g in e['region'])),flush=True)
def panels(f,entry,rects):
 img=entry['img'];xy=entry['xy'];expr=entry['expr'];region=entry['region'];assert len(xy)==len(expr)==len(region) and np.isfinite(expr).all()
 for j,rect in enumerate(rects):
  ax=f.add_axes(rect);ax.imshow(img,alpha=1 if j==0 else .52);ax.set_aspect('equal');ax.set_xlim(0,img.shape[1]);ax.set_ylim(img.shape[0],0);ax.set_xticks([]);ax.set_yticks([])
  for spine in ax.spines.values():spine.set_visible(False)
  if j==1:
   for group in sorted(set(region)):
    m=region==group;ax.scatter(xy[m,0],xy[m,1],s=3.2,c=COL[group],lw=0,alpha=.88,rasterized=True)
  if j==2:
   zero=expr==0;ax.scatter(xy[zero,0],xy[zero,1],s=3.2,c=ZERO,lw=0,alpha=.88,rasterized=True)
   order=np.argsort(expr,kind='stable');pos=order[expr[order]>0];ax.scatter(xy[pos,0],xy[pos,1],s=3.2,c=expr[pos],cmap='viridis',vmin=0,vmax=CAP,lw=0,alpha=1,rasterized=True)
def save(f,name):
 for ext in ['png','pdf','svg']:f.savefig(P/(name+'.'+ext),dpi=250,facecolor='white')
 plt.close(f)
for e in data:
 f=plt.figure(figsize=(12,4.4));rects=[[x,.16,.27,.71] for x in [.025,.328,.631]];panels(f,e,rects)
 for x,t in zip([.025,.328,.631],['H&E','作者病理区域','LYPLA1 实测表达']):f.text(x,.92,t,fontsize=14,color=G,weight='bold')
 cats=[g for g in COL if g in set(e['region'])];handles=[Line2D([0],[0],marker='o',ls='',color=COL[g],label=LABEL[g],markersize=5) for g in cats]
 f.legend(handles=handles,loc='lower center',bbox_to_anchor=(.48,.02),ncol=len(cats),frameon=False,fontsize=9)
 cb=f.add_axes([.936,.28,.011,.46]);f.colorbar(ScalarMappable(norm=Normalize(0,CAP),cmap='viridis'),cax=cb,extend='max',ticks=[0,1,2,3]);cb.tick_params(labelsize=9);cb.set_ylabel('log1p(CP10K)',fontsize=9)
 f.suptitle(e['cancer']+' · '+e['sample']+' · '+e['patient'],fontsize=12,y=.995,color=M)
 save(f,e['cancer']+'_'+e['sample'])
 checks.append(dict(cancer=e['cancer'],sample=e['sample'],patient=e['patient'],n_spots=len(e['expr']),n_positive=int((e['expr']>0).sum()),regions={LABEL[g]:int((e['region']==g).sum()) for g in set(e['region'])},cap=CAP,n_above_cap=int((e['expr']>CAP).sum())))
for page,c in [(15,'BRCA'),(16,'PRAD')]:
 f=plt.figure(figsize=(16,9));entries=[e for e in data if e['cancer']==c]
 def text(x,y,s,size=12,col=INK,bold=False,ha='left'):f.text(x,y,s,fontsize=size,color=col,weight='bold' if bold else 'normal',ha=ha,va='center')
 def line(y):f.add_artist(Line2D([.05,.95],[y,y],transform=f.transFigure,color=LINE,lw=.8))
 ax=f.add_axes([.044,.882,.168,.10]);ax.imshow(plt.imread(AS/'sysu_logo.png'));ax.axis('off');text(.95,.93,'CAMP  /  组会汇报',10,M,ha='right');line(.883)
 title='BRCA：LYPLA1 在病理癌区中的空间分布' if c=='BRCA' else 'PRAD：对照癌区与良性腺体的 LYPLA1 表达'
 text(.05,.828,title,26,G,True)
 sub='CTA · 两位患者各一张代表切片 · 癌区与免疫区采用作者病理标注' if c=='BRCA' else 'Erickson 2022 · 同一患者的两张代表切片 · 作者病理共识标注'
 text(.05,.776,sub,11.5,M)
 for x,t in zip([.16,.42,.68],['H&E','病理区域','LYPLA1 实测表达']):text(x,.720,t,14,G,True)
 for e,y in zip(entries,[.420,.093]):
  panels(f,e,[[x,y,.225,.286] for x in [.16,.42,.68]])
  sample=e['sample'].replace('V19T26-','')
  text(.05,y+.259,sample,12,G,True);text(.05,y+.226,e['patient'],10,M)
  cats=[g for g in COL if g in set(e['region'])]
  for k,g in enumerate(cats):
   yy=y+.179-k*.027;f.add_artist(Line2D([.057],[yy],transform=f.transFigure,marker='o',ls='',color=COL[g],markersize=4.8));text(.068,yy,LABEL[g],9,INK)
 line(.402)
 cb=f.add_axes([.925,.285,.008,.28]);f.colorbar(ScalarMappable(norm=Normalize(0,CAP),cmap='viridis'),cax=cb,extend='max',ticks=[0,1,2,3]);cb.tick_params(labelsize=9,length=2);cb.set_ylabel('LYPLA1 · log1p(CP10K)',fontsize=9,color=M,labelpad=6)
 line(.076)
 foot='未标注区不等于正常上皮；代表切片用于定位，不替代完整患者汇总。' if c=='BRCA' else '两张切片均来自 Patient 1，不代表两位患者；多细胞 spot 不等于纯癌细胞。'
 text(.05,.050,foot,10,G,True);text(.05,.024,'统一色标 0–3；浅灰为零表达，超出上限按端点色显示。不同研究不直接比较绝对表达。',8.3,M);text(.96,.030,str(page),13,G,True,ha='right')
 save(f,f'{page}_{c}_LYPLA1_空间切片统一版')
(P/'source_manifest.json').write_text(json.dumps(sources,indent=2))
(P/'display_checks.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2))
(P/'analysis_spec.json').write_text(json.dumps(dict(selected_by='user screenshots',new_statistics=False,source_coordinates_reused=True,source_annotations_reused=True,cap=CAP,cmap='viridis',zero_color=ZERO,marker_area=3.2,HE_alpha=1,overlay_HE_alpha=.52,scale='log1p(CP10K)',interpolation=False),indent=2))
print('DONE: four individual figures and two slides',flush=True)
