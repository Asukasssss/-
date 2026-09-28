"""All687 descriptive dotplots and heatmaps from public donor-equal profiles only."""
import json,math,html
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from record_sop_v3_resume import SOURCE,MAP,write,checksum,read

def main():
    df=read(SOURCE/'sc_celltype_profiles.tsv');pool=pd.read_csv(MAP/'gene_pool_history_union.tsv',sep='\t');genes=sorted(pool.gene);cohorts=['GSE263733','GSE278688','GSE242230'];out=SOURCE/'figures';out.mkdir(exist_ok=True)
    assert len(genes)==687 and set(df.gene)==set(genes) and not df.duplicated(['cohort','gene','celltype']).any()
    globalmax=float(df.effect.max());norm=Normalize(0,globalmax);cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#d5d9de')
    identity_v4=(SOURCE/'author_identity_counts.tsv').exists()
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    manifest=[];npages=math.ceil(len(genes)/48)
    for page in range(npages):
        gg=genes[page*48:(page+1)*48]
        for kind in ['dotplot','heatmap']:
            fig,axes=plt.subplots(1,3,figsize=(18,15),sharey=True,gridspec_kw={'wspace':.18})
            fig.subplots_adjust(left=.085,right=.92,top=.89,bottom=.21 if identity_v4 else .17)
            for ax,cohort in zip(axes,cohorts):
                d=df[df.cohort==cohort];cats=sorted(d.celltype.unique());x=d.pivot(index='gene',columns='celltype',values='effect').reindex(index=gg,columns=cats).to_numpy(float);det=d.pivot(index='gene',columns='celltype',values='mean_detection_fraction').reindex(index=gg,columns=cats).to_numpy(float)
                if kind=='heatmap':ax.imshow(np.ma.masked_invalid(x),aspect='auto',interpolation='nearest',cmap=cmap,norm=norm)
                else:
                    ax.set_facecolor('#fafafa');yy,xx=np.indices(x.shape);ok=np.isfinite(x)
                    ax.scatter(xx[~ok],yy[~ok],s=25,marker='s',color='#d5d9de',edgecolors='none')
                    ax.scatter(xx[ok],yy[ok],s=4+110*det[ok],c=x[ok],cmap=cmap,norm=norm,edgecolors='#6b7280',linewidths=.15)
                    ax.set_ylim(len(gg)-.5,-.5);ax.set_xlim(-.5,len(cats)-.5)
                    ax.set_xticks(np.arange(len(cats))-.5,minor=True);ax.set_yticks(np.arange(len(gg))-.5,minor=True);ax.grid(which='minor',color='#e6e9ec',linewidth=.3);ax.tick_params(which='minor',length=0)
                ax.set_xticks(np.arange(len(cats)));ax.set_xticklabels(cats,rotation=65,ha='right',fontsize=8);ax.set_yticks(np.arange(len(gg)));ax.set_yticklabels(gg,fontsize=8);ax.tick_params(axis='y',length=0);ax.set_title(cohort,fontsize=12,pad=9)
            barax=fig.add_axes([.94,.43,.012,.32]);fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=barax,label='Equal-label mean cellwise log1p(counts per 10k)')
            fig.suptitle(f'PDAC | All-candidate {kind} | {page+1}/{npages}',x=.085,ha='left',y=.965,fontsize=18,weight='bold')
            fig.text(.085,.931,f'{gg[0]} - {gg[-1]} | Alphabetical pages; all687 genes retained. No RNA/association significance filter.',fontsize=10)
            fig.text(.085,.03,('Author identities: malignant, normal epithelium, unresolved ductal are distinct. ' if identity_v4 else 'Author broad cell labels. ')+'>=20cells per label/type;>=3labels per type. Gray = missing, not zero.\nSource labels equally weighted. GSE242230 uses author sample labels. Descriptive expression only; no SC P/q or mechanism inference.',fontsize=9,linespacing=1.5)
            if kind=='dotplot':
                hs=[axes[2].scatter([],[],s=4+110*v,facecolor='#718096',edgecolor='#718096') for v in [0,.1,.5,1]]
                fig.legend(hs,['0%','10%','50%','100%'],title='Equal-label detection fraction (dot area)',loc='lower center',bbox_to_anchor=(.59,.053),ncol=4,frameon=False,fontsize=9,title_fontsize=9)
            stem=f'{kind}_{page+1:02d}';fig.savefig(out/(stem+'.svg'));fig.savefig(out/(stem+'.png'),dpi=110);plt.close(fig)
            svg=out/(stem+'.svg');svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n',encoding='utf8',newline='\n')
            manifest.append({'file':stem+'.svg','preview':stem+'.png','kind':kind,'page':page+1,'genes':';'.join(gg),'n_genes':len(gg),'source_table':'../sc_celltype_profiles.tsv','selection':'All687 alphabetical;48genes/page','color_min':0,'color_max':globalmax,'detection_area':'4+110*fraction for dotplots;zero retains tiny marker','meaning':'Donor-equal expression source;no disease comparison or mechanism'})
    write(out/'figure_manifest.tsv',pd.DataFrame(manifest))
    text=['# PDAC 全候选来源图','',f'覆盖687基因，每页最多48个，共{npages}页点图及{npages}页热图。全部按字母排序，不按P筛选。','', '颜色为先逐细胞全库log1p(每万计数)，再来源标签内均值、标签等权均值。点面积表示标签等权检出比例（保留小的零值标记）。灰色为未测或覆盖不足，不表示0。所有页使用同一颜色范围；不同研究绝对表达差异不用于正式研究间检验。GSE242230单位是作者样本标签。','', '所有原矩阵、逐细胞坐标及逐供者数值留服务器。图只使用公开的供者等权汇总；无单细胞P/q、无机制推断。','']
    for i in range(npages):text += [f'## 第{i+1}页','',f'![点图](dotplot_{i+1:02d}.png)','',f'[点图SVG](dotplot_{i+1:02d}.svg) · [热图PNG](heatmap_{i+1:02d}.png) · [热图SVG](heatmap_{i+1:02d}.svg)','']
    (out/'README_CN.md').write_text('\n'.join(text),encoding='utf8',newline='\n')
    # Compact browser index; original SVG remains independently usable.
    cards=''.join(f'<article><h2>{kind} {i+1}</h2><a href="{kind}_{i+1:02d}.svg"><img loading="lazy" src="{kind}_{i+1:02d}.png"></a></article>' for i in range(npages) for kind in ['dotplot','heatmap'])
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>PDAC expression sources</title><style>body{font:16px system-ui;margin:30px;background:#f3f5f7}main{display:grid;grid-template-columns:1fr 1fr;gap:20px}article{background:white;padding:15px}img{width:100%}</style><h1>PDAC全687基因表达来源</h1><p>全部候选按字母分15页。灰色为不可评估；图形仅描述来源，不证明机制。点击查看矢量图。</p><main>'+cards+'</main>',encoding='utf8')
    write_json={'status':'DONE','genes':687,'pages_per_kind':npages,'figures':len(manifest),'normalization':'donor equal mean of cellwise log1p10k','no_significance_selection':True,'matplotlib':matplotlib.__version__}
    (out/'plot_spec.json').write_text(json.dumps(write_json,indent=2)+'\n');print(json.dumps(write_json))

if __name__=='__main__':main()
