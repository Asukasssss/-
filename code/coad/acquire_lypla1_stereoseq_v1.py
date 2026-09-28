"""Download public STT0000036 data to server165 only; metadata precedes expression."""
import argparse,concurrent.futures,hashlib,json,urllib.request,time
from pathlib import Path
import pandas as pd
DATA=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/coad_lypla1_spatial_alternative_20260928')
def get(u):
 p=DATA/u.split('/')[-1]
 if not p.exists():
  for k in range(3):
   try:
    with urllib.request.urlopen(u,timeout=90)as r,open(str(p)+'.part','wb')as f:
     while True:
      b=r.read(1024*1024)
      if not b:break
      f.write(b)
    Path(str(p)+'.part').replace(p);break
   except Exception:
    if k==2:raise
    time.sleep(3)
 h=hashlib.sha256()
 with open(p,'rb')as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return dict(url=u,file=p.name,bytes=p.stat().st_size,sha256=h.hexdigest())
def main():
 a=argparse.ArgumentParser();a.add_argument('catalog');a.add_argument('--expression',action='store_true');v=a.parse_args()
 rows=json.load(open(v.catalog));urls=[u for r in rows for u in r['files'] if u.endswith('.csv')]
 with concurrent.futures.ThreadPoolExecutor(3)as ex:proof=list(ex.map(get,urls))
 overview=[]
 for u in urls:
  d=pd.read_csv(DATA/u.split('/')[-1]);x=dict(file=u.split('/')[-1],n=len(d))
  for k in ['id','Type','Treatment','group1','group2']:
   assert d[k].nunique(dropna=False)==1;x[k]=str(d[k].iloc[0])
  overview.append(x)
 pd.DataFrame(overview).to_csv(DATA/'section_catalog.tsv',sep='\t',index=False)
 print(json.dumps(overview),flush=True)
 if v.expression:
  admitted={r['id'] for r in overview if r['Treatment'].lower() in ('none','no','untreated','naive','nan','na')}
  assert admitted,'Inspect actual treatment labels first'
  urls=[u for r in rows for u in r['files'] if u.endswith('.gem.gz') and u.split('/')[-1][13:-7] in admitted]
  assert len(urls)==len(admitted)
  with concurrent.futures.ThreadPoolExecutor(3)as ex:
   for p in ex.map(get,urls):proof.append(p);print('READY',p['file'],p['bytes'],flush=True)
 json.dump(proof,open(DATA/'stereoseq_download_manifest.json','w'),indent=2)
if __name__=='__main__':main()
