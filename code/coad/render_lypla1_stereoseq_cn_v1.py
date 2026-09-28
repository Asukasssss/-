"""Chinese reading guide; preserves expression and coordinates, adds no statistics."""
import json
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
RUN=Path(__file__).resolve().parent;OUT=RUN/'public'
font_manager.fontManager.addfont('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
plt.rcParams['font.family']='Noto Sans CJK JP';plt.rcParams['axes.unicode_minus']=False
NAMES={'Tumor-associated epithelium':'癌相关上皮区（作者注释）','Normal epithelium':'正常上皮区（作者注释）','Tumor-stroma boundary':'肿瘤与基质交界区','Stroma':'基质区','Immune-rich stroma':'免疫细胞丰富区','Muscle':'肌层','Fiber/cavity':'纤维／空腔区域','Other':'其他'}
COLORS={'Tumor-associated epithelium':'#d6334a','Normal epithelium':'#2379bd','Tumor-stroma boundary':'#e3af28','Stroma':'#9ba7b4','Immune-rich stroma':'#49a582','Muscle':'#9671ad','Fiber/cavity':'#e2e2e2','Other':'#555555'}
spec=json.load(open(OUT/'analysis_spec.json'));vmax=spec['color_range'][1]
cmap=LinearSegmentedColormap.from_list('lypla1',['#eeeeee','#ffcc80','#ed7038','#a60026'])
for p in sorted(RUN.glob('*_private.tsv.gz')):
 sid=p.name.replace('_private.tsv.gz','');m=pd.read_csv(p,sep='\t')
 fig,axs=plt.subplots(1,2,figsize=(12,7.6));fig.suptitle(sid+'｜LYPLA1空间表达',fontsize=19,y=.97)
 fig.text(.5,.92,'同一张切片、相同坐标：先看左边的组织标签，再看右边对应位置的表达',ha='center',fontsize=11)
 for reg,c in COLORS.items():
  d=m[m.region.eq(reg)];axs[0].scatter(d.x,d.y,c=c,s=5,marker='s',linewidths=0,rasterized=True)
 z=axs[1].scatter(m.x,m.y,c=m.LYPLA1,cmap=cmap,vmin=0,vmax=vmax,s=5,marker='s',linewidths=0,rasterized=True)
 xs=np.arange(m.x.min()-1,m.x.max()+2);ys=np.arange(m.y.min()-1,m.y.max()+2);mask=np.zeros((len(ys),len(xs)))
 t=m[m.region.eq('Tumor-associated epithelium')];mask[(t.y-ys[0]).astype(int),(t.x-xs[0]).astype(int)]=1
 if len(t):axs[1].contour(xs,ys,mask,levels=[.5],colors=['#111111'],linewidths=.4)
 axs[0].set_title('① 红色在哪里：癌相关上皮区',fontsize=13);axs[1].set_title('② 颜色越深：LYPLA1表达越高',fontsize=13)
 for ax in axs:ax.set_aspect('equal');ax.invert_yaxis();ax.axis('off')
 cb=fig.colorbar(z,ax=axs[1],fraction=.035,pad=.02);cb.set_label('LYPLA1：log1p(CP10K)',fontsize=10)
 present=[r for r in COLORS if r in set(m.region)]
 fig.legend(handles=[Patch(color=COLORS[r],label=NAMES[r])for r in present],loc='lower center',ncol=3,fontsize=10,bbox_to_anchor=(.5,.09))
 a=m[m.region.eq('Tumor-associated epithelium')].LYPLA1;b=m[m.region.eq('Normal epithelium')].LYPLA1
 txt='';
 if len(a) and len(b):txt='本切片平均表达：癌相关上皮 %.3f；正常上皮 %.3f（描述性比较）'%(a.mean(),b.mean())
 fig.text(.5,.062,txt,ha='center',fontsize=11)
 fig.text(.5,.025,'右图黑线对应左图癌相关上皮区；红色表达不等于癌区。每个方点为约50微米空间单元，可含多种细胞。',ha='center',fontsize=9)
 fig.subplots_adjust(left=.02,right=.96,bottom=.27,top=.84,wspace=.13)
 for ext in ['png','pdf']:fig.savefig(OUT/'figures'/(sid+'_LYPLA1_CN.'+ext),dpi=180,bbox_inches='tight',pad_inches=.15)
 plt.close(fig)
print('CHINESE_FIGURES_DONE',flush=True)
