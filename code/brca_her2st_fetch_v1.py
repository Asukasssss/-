"""Pinned author repository files, RAM-only relay to server165."""
import requests,base64,hashlib,json,concurrent.futures,time
from brca_gse210616_stream_to_server import transfer,remote
BASE='/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/BRCA_Andersson_HER2ST'
API='https://api.github.com/repos/almaan/her2st'
def get(url):
 for i in range(4):
  try:
   r=requests.get(url,timeout=45);r.raise_for_status();return r.json()
  except Exception:
   if i==3:raise
   time.sleep(2)
tree=get(API+'/git/trees/master?recursive=1');sha=tree['sha']
samples=['A1','B1','C1','D1','E1','F1','G2','H1']
paths=set(['README.md','data/ST-pat/lbl/README'])
for s in samples:
 paths.update([f'data/ST-cnts/{s}.tsv.gz',f'data/ST-pat/lbl/{s}_labeled_coordinates.tsv',f'data/ST-spotfiles/{s}_selection.tsv'])
items=[x for x in tree['tree'] if x['type']=='blob' and (x['path'] in paths or any(x['path'].startswith(f'data/ST-imgs/{s[0]}/{s}/') for s in samples))]
def one(x):
 p=BASE+'/'+x['path'];remote("mkdir -p '"+p.rsplit('/',1)[0]+"'")
 b=base64.b64decode(get(API+'/git/blobs/'+x['sha'])['content'])
 assert hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==x['sha']
 transfer([b],p);print(x['path'],len(b),flush=True)
 return {'path':p,'url':f'https://github.com/almaan/her2st/blob/{sha}/'+x['path'],'git_blob_sha':x['sha'],'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b)}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:records=list(ex.map(one,items))
transfer([json.dumps({'commit_tree':sha,'files':records},indent=2).encode()],BASE+'/manifest.json')
print('DOWNLOAD_DONE',len(records),flush=True)
