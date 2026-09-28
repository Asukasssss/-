"""Register author pathology SVG images to Space Ranger images, blinded to LYPLA1."""
import pathlib,xml.etree.ElementTree as ET,base64,io,re,json,collections
import numpy as np
from PIL import Image
from skimage.color import rgb2gray
from skimage.feature import SIFT,match_descriptors
from skimage.measure import ransac
from skimage.transform import AffineTransform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
D=pathlib.Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/BRCA_CTA2025')
O=pathlib.Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A/20260928T060000Z_cta_lypla1_v1')
XL='{http://www.w3.org/1999/xlink}'
def transform(s):
 if not s:return np.eye(3)
 m=re.fullmatch(r'matrix\(([^)]+)\)',s)
 if not m:raise ValueError(s)
 a,b,c,d,e,f=map(float,re.split('[ ,]+',m[1]));return np.array([[a,c,e],[b,d,f],[0,0,1.]])
def coords(s):
 toks=re.findall(r'[a-zA-Z]|[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?',s)
 i=0;cmd=None;cur=np.zeros(2);first=None;out=[];allpaths=[]
 sizes={'M':2,'L':2,'H':1,'V':1,'C':6,'Q':4,'Z':0}
 while i<len(toks):
  if toks[i].isalpha():cmd=toks[i];i+=1
  up=cmd.upper()
  if up not in sizes:raise ValueError('Unsupported SVG command '+cmd)
  if up=='Z':
   if len(out)>2:allpaths.append(np.array(out))
   cur=first.copy();out=[];cmd=None;continue
  n=sizes[up];v=np.array(list(map(float,toks[i:i+n])));i+=n;rel=cmd.islower()
  if up in ['M','L']:
   cur=v+(cur if rel else 0)
   if up=='M':first=cur.copy();cmd='l' if rel else 'L'
   out.append(cur.copy())
  elif up in ['H','V']:
   k=0 if up=='H' else 1;cur[k]=v[0]+(cur[k] if rel else 0);out.append(cur.copy())
  else:
   ps=v.reshape(-1,2)+(cur if rel else 0);t=np.linspace(0,1,9)[1:,None]
   xy=(1-t)**3*cur+3*(1-t)**2*t*ps[0]+3*(1-t)*t*t*ps[1]+t**3*ps[2] if up=='C' else (1-t)**2*cur+2*(1-t)*t*ps[0]+t*t*ps[1]
   out.extend(xy);cur=ps[-1].copy()
 return allpaths
def parse(p):
 root=ET.parse(p).getroot();images={e.attrib['id']:e for e in root.iter() if e.tag.endswith('image')};uses=[];polys=[];styles=collections.Counter()
 def visit(e,parent):
  mat=parent@transform(e.attrib.get('transform',''))
  if e.tag.endswith('use') and e.attrib.get(XL+'href','')[1:] in images:uses.append((e,mat,images[e.attrib[XL+'href'][1:]]))
  if e.tag.endswith('path') and 'd' in e.attrib:
   style=e.attrib.get('style','');m=re.search(r'stroke:(#[a-fA-F0-9]+)',style)
   if m:
    color=m[1].lower();styles[color]+=1
    h=color[1:];h=''.join(x*2 for x in h) if len(h)==3 else h
    rgb=np.array([int(h[j:j+2],16) for j in [0,2,4]])
    label='Tumor' if rgb[0]>150 and rgb[1]<80 and rgb[2]<80 else 'Immune' if rgb[0]>180 and rgb[1]>120 and rgb[2]<90 else 'DCIS' if rgb[2]>100 and rgb[0]<80 and rgb[1]<100 else 'Vessel' if rgb[1]>120 and rgb[2]>120 and rgb[0]<150 else 'Necrosis' if max(rgb)<100 else 'Unknown'
    for xy in coords(e.attrib['d']):polys.append({'label':label,'color':color,'xy':(np.c_[xy,np.ones(len(xy))]@mat.T)[:,:2]})
  for c in e:visit(c,mat)
 visit(root,np.eye(3));assert len(uses)==1
 use,mat,im=uses[0];pic=Image.open(io.BytesIO(base64.b64decode(im.attrib[XL+'href'].split(',',1)[1]))).convert('RGB')
 # SVG image intrinsic dimensions and the <use> viewport coincide in this source.
 assert abs(float(re.sub('px','',use.attrib['width']))-pic.width)<1 and abs(float(re.sub('px','',use.attrib['height']))-pic.height)<1
 return np.array(pic),mat,polys,dict(styles)
def features(im):
 f=SIFT(upsampling=1);f.detect_and_extract(rgb2gray(im));return f.keypoints[:,::-1],f.descriptors
def main():
 O.mkdir(exist_ok=True,parents=True);(O/'registration').mkdir(exist_ok=True)
 targets={}
 for p in D.glob('spaceranger_output/*/outs/spatial/tissue_lowres_image.png'):
  im=np.array(Image.open(p).convert('RGB'));targets[p.parents[2].name]=(im,*features(im))
 results=[]
 for p in sorted(D.glob('Images/Manual_annotation/*.svg')):
  try:
   im,mat,polys,styles=parse(p);kp,desc=features(im);hits=[]
   # Restrict by explicitly matching slide barcode, then register sections by image content.
   card=re.search(r'V\d\d[A-Z]\d\d-\d\d\d',p.name)
   explicit_area=re.search(r'\.([A-D]1)(?:_|\.)',p.name)
   for name,(dest,qk,qd) in targets.items():
    if card and not name.startswith(card[0]):continue
    if explicit_area and name != card[0]+'_'+explicit_area[1]:continue
    matches=match_descriptors(desc,qd,cross_check=True,max_ratio=.7)
    if len(matches)<8:continue
    src=kp[matches[:,0]];dst=qk[matches[:,1]]
    model,inside=ransac((src,dst),AffineTransform,min_samples=3,residual_threshold=2,max_trials=2000,rng=42)
    if model is None:continue
    nin=int(inside.sum());err=float(np.median(np.linalg.norm(model(src[inside])-dst[inside],axis=1)))
    hits.append((nin,name,model,err,len(matches)))
   hits.sort(key=lambda x:x[0],reverse=True)
   if not hits:results.append({'svg':p.name,'status':'NOT_EVALUABLE','reason':'No registration match'});continue
   nin,name,model,err,nm=hits[0];second=hits[1][0] if len(hits)>1 else 0
   passed=nin>=20 and nin/nm>=.35 and err<=1.5 and second<nin*.65
   t=model.params@np.linalg.inv(mat)
   aligned=[dict(label=x['label'],color=x['color'],xy=(np.c_[x['xy'],np.ones(len(x['xy']))]@t.T)[:,:2].tolist()) for x in polys]
   row={'svg':p.name,'sample':name,'n_inliers':nin,'n_matches':nm,'median_error_lowres_px':err,'second_match_inliers':second,'status':'PASS' if passed else 'NEEDS_REVIEW','transform_svg_to_lowres':t.tolist(),'styles':styles}
   results.append(row)
   (O/'registration'/f'{name}__{p.stem}.json').write_text(json.dumps({'registration':row,'polygons':aligned}))
   fig,ax=plt.subplots(figsize=(8,7));ax.imshow(targets[name][0])
   for x in aligned:
    xy=np.array(x['xy']);ax.plot(*np.r_[xy,xy[:1]].T,color=x['color'],lw=.6)
   ax.set_title(f'{name} | {nin} inliers | median error {err:.2f}px');ax.set_xlim(0,targets[name][0].shape[1]);ax.set_ylim(targets[name][0].shape[0],0);ax.axis('off');fig.savefig(O/'registration'/f'{name}__{p.stem}.png',dpi=140,bbox_inches='tight');plt.close(fig)
   print(p.name,name,row['status'],nin,err,flush=True)
  except Exception as e:
   results.append({'svg':p.name,'status':'NEEDS_REVIEW','reason':repr(e)});print(p.name,repr(e),flush=True)
 (O/'registration_summary.json').write_text(json.dumps(results,indent=2))
if __name__=='__main__':main()
