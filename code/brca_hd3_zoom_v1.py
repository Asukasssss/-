"""Render a tissue-focused view from saved 16um measurements; no statistical rerun."""
from pathlib import Path
import json,sys,io,base64
import numpy as np,pandas as pd
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(sys.argv[1]);R=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/BRCA_HD3_Kinnex')
d=pd.read_csv(O/'bin_16um_SERVER_ONLY.tsv.gz',sep='\t');j=json.JSONDecoder().raw_decode((R/'web_summary.html').read_text().split('const data = ',1)[1])[0]
p=j['tabs']['tab_data'][2]['clustering']['card']['inner']['per_bin_level']['options'][0]['component']['plot']['spatial_plot_props']
im=np.array(Image.open(io.BytesIO(base64.b64decode(j['_resources']['000'].split(',',1)[1]))))
T=np.linalg.inv(np.array(p['tissue_css_transform']).reshape(4,4,order='F'))@np.array(p['spot_css_transform']).reshape(4,4,order='F')
v=np.vstack([(d.bin_col.to_numpy()+.5)*2,(d.bin_row.to_numpy()+.5)*2,np.zeros(len(d)),np.ones(len(d))]);xy=T@v;xy=xy[:2]/xy[3]
good=d.qc.to_numpy();total=d.total_counts.to_numpy();target=np.log1p(d.LYPLA1_count.to_numpy()/np.maximum(total,1)*1e4)
epi=np.mean([np.log1p(d[g+'_count'].to_numpy()/np.maximum(total,1)*1e4) for g in ['EPCAM','KRT8','KRT18','KRT19']],axis=0)
x0,x1=np.quantile(xy[0,good],[0,1]);y0,y1=np.quantile(xy[1,good],[0,1]);pad=.07*max(x1-x0,y1-y0)
fig,axs=plt.subplots(1,3,figsize=(12,7))
for ax in axs:ax.imshow(im);ax.set_xlim(x0-pad,x1+pad);ax.set_ylim(y1+pad,y0-pad);ax.axis('off')
axs[0].set_title('H&E; no pathological region labels')
for ax,title,y in [(axs[1],'LYPLA1',target),(axs[2],'Epithelial marker mean',epi)]:
 s=ax.scatter(xy[0,good],xy[1,good],s=.9,c=y[good],cmap='magma',vmin=0,vmax=max(.1,np.quantile(y[good],.995)),alpha=.85,rasterized=True);ax.set_title(title);fig.colorbar(s,ax=ax,shrink=.55,label='log1p(CP10K)')
fig.suptitle('Breast IDC | Visium HD 3-prime + Kinnex | 16 um\nOne specimen; epithelial marker signal does not identify malignancy')
fig.tight_layout();fig.savefig(O/'figures'/'he_zoom_16um.png',dpi=200);plt.close(fig)
