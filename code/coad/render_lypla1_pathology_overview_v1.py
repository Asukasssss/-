"""Layout-only overview: reuse frozen spot measurements; no statistical recomputation."""
import argparse,hashlib,io,json,zipfile
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.collections import PatchCollection
from matplotlib.patches import Patch,Circle
from PIL import Image
from acquire_lypla1_pathology_v1 import DATA,SAMPLES
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
pub=a.out/'public';fig=pub/'figures'
v=json.loads((pub/'validation.json').read_text());norm=Normalize(*v['shared_expression_range'])
cmap=LinearSegmentedColormap.from_list('LYPLA1',['#f0f0f0','#ffd074','#ef703c','#980019'])
colors={'Tumor':'#cf243c','Tumor-stroma mixed':'#f39935','Stroma':'#237ec5','Non-neoplastic epithelium':'#2b9c62','Other annotated tissue':'#a481b1','Excluded / unannotated / QC fail':'#bdbdbd'}
ff,axs=plt.subplots(7,2,figsize=(10,25))
for i,(alias,name,_) in enumerate(SAMPLES):
 f=pd.read_csv(a.out/'private'/(alias+'_spots.tsv'),sep='\t');f=f[f.in_tissue==1]
 with zipfile.ZipFile(DATA/(name+'.selected.zip')) as z:
  def get(s):return z.read(next(n for n in z.namelist() if n.endswith(s)))
  img=np.asarray(Image.open(io.BytesIO(get('tissue_lowres_image.png'))));scale=json.loads(get('scalefactors_json.json'))
 r=scale['spot_diameter_fullres']*scale['tissue_lowres_scalef']/2
 margin=15;lims=(max(0,f.x.min()-margin),min(img.shape[1],f.x.max()+margin),max(0,f.y.min()-margin),min(img.shape[0],f.y.max()+margin))
 for j,ax in enumerate(axs[i]):
  ax.imshow(img,alpha=.38)
  col=[colors[s] for s in f.display_group] if j==0 else [cmap(norm(n)) if ok else '#bdbdbd' for n,ok in zip(f.log1pCP10K,f.admitted)]
  ax.add_collection(PatchCollection([Circle((x,y),r*.92) for x,y in zip(f.x,f.y)],facecolors=col,edgecolors='none',alpha=.95,rasterized=True))
  if j==1:
   tumor=f[f.admitted&(f.original_broad=='Tumor')]
   ax.add_collection(PatchCollection([Circle((x,y),r*1.04) for x,y in zip(tumor.x,tumor.y)],facecolors='none',edgecolors='#00a4b5',linewidths=.55,rasterized=True))
  ax.set_xlim(lims[0],lims[1]);ax.set_ylim(lims[3],lims[2]);ax.set_aspect('equal');ax.axis('off')
  ax.set_title(alias+' | pathology' if j==0 else 'LYPLA1 | cyan rings = tumor',fontsize=11)
ff.legend(handles=[Patch(color=c,label=l) for l,c in colors.items()],loc='lower center',ncol=2,fontsize=10)
ff.subplots_adjust(bottom=.055,top=.945,left=.04,right=.88,wspace=.12,hspace=.22)
cb=ff.add_axes([.92,.4,.018,.2]);ff.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cb,label='LYPLA1 log1p(CP10K)')
ff.suptitle('CRC: pathologist regions and LYPLA1',fontsize=16,y=.985)
ff.text(.5,.97,'Matched coordinates; shared expression scale; one section per case',ha='center',fontsize=10)
for ext in ['png','pdf']:ff.savefig(fig/('all7_pathology_LYPLA1.'+ext),dpi=150,bbox_inches='tight')
plt.close(ff)
v['overview_render_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
v['overview_layout']='Tissue-area view limits; original coordinates and all in-tissue dots preserved; titles separated; no numeric change'
(pub/'validation.json').write_text(json.dumps(v,indent=2)+'\n')
print('OVERVIEW_RENDERED')
