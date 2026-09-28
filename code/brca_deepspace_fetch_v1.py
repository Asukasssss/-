"""Retrieve public expert polygons, preserving source data on server165 only."""
import requests,json,time,hashlib,concurrent.futures
from brca_gse210616_stream_to_server import transfer,remote
B='https://genomics.virus.kyoto-u.ac.jp/deepspacedb'
D='/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/BRCA_DeepSpaceDB'
def get(url,**kwargs):
 for n in range(3):
  try:
   r=requests.get(url,timeout=45,**kwargs)
   if r.status_code==429: raise RuntimeError('Service rate limit: stop requests and use original public archive')
   r.raise_for_status();return r
  except RuntimeError: raise
  except Exception:
   if n==2:raise
   time.sleep(2)
def one(x):
 u=B+'/database/'+x['ID']+'/annotation_manual_plot';r=get(u);j=r.json()
 p=json.loads(j['manual_annotation_plot']) if j.get('manual_annotation_plot') else None
 names=sorted(set(t['name'] for t in p['data'])) if p else []
 print(x['ID'],x['series_id'],x['sample_id'],names,flush=True)
 return {'metadata':x,'url':u,'sha256':hashlib.sha256(r.content).hexdigest(),'response':j,'labels':names}
if __name__=='__main__':
 remote("mkdir -p '"+D+"'")
 r=get(B+'/get_data',params={'search':'breast','length':500,'start':0,'draw':1});j=r.json();assert len(j['data'])==j['recordsFiltered']
 selected=[x for x in j['data'] if x.get('has_annotation')=='yes' and x.get('organism')=='human']
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex: rec=list(ex.map(one,selected))
 transfer([json.dumps({'metadata_url':r.url,'records':rec},ensure_ascii=False).encode()],D+'/expert_annotations.json')
