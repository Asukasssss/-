"""GSE213688 pathology-first eligibility audit, RAM relay of source files to server165."""
import requests,re,json,concurrent.futures,time,hashlib,sys
from brca_gse210616_stream_to_server import transfer,remote
BASE='/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/BRCA_GSE213688'
def get(u):
 for a in range(4):
  try:r=requests.get(u,timeout=45);r.raise_for_status();return r
  except Exception:
   if a==3:raise
   time.sleep(2)
def send(b,p):
 for a in range(4):
  try:transfer([b],p);return
  except Exception:
   if a==3:raise
   time.sleep(2)
mode=sys.argv[1] if len(sys.argv)>1 else 'annotations'
remote("mkdir -p '"+BASE+"'")
u='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213688&targ=self&form=text&view=quick';r=get(u);send(r.content,BASE+'/series.soft')
gsms=re.findall(r'!Series_sample_id = (GSM\d+)',r.text);assert len(gsms)==15
def one(gsm):
 d=BASE+'/'+gsm;remote("mkdir -p '"+d+"'");r=get(f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=quick');send(r.content,d+'/metadata.soft')
 urls=re.findall(r'!Sample_supplementary_file(?:_\d+)? = (\S+)',r.text);title=re.search(r'!Sample_title = (.+)',r.text).group(1).strip()
 rec={'gsm':gsm,'title':title,'characteristics':re.findall(r'!Sample_characteristics_ch1 = (.+)',r.text),'files':[]}
 for url in urls:
  if mode=='annotations' and 'pathologist_annotations' not in url:continue
  if mode=='all' and not any(k in url for k in ['pathologist_annotations','filtered_feature_bc_matrix','scalefactors','tissue_hires_image','tissue_positions']):continue
  url=url.replace('ftp://','https://');data=get(url).content;p=d+'/'+url.rsplit('/',1)[1];send(data,p)
  rec['files'].append({'path':p,'url':url,'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data)})
 print(gsm,title,'DONE',flush=True);return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:records=list(ex.map(one,gsms))
send(json.dumps(records,indent=2).encode(),BASE+'/manifest_'+mode+'.json')
