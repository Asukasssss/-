"""Fetch author's public scRNA data directly to server165; no local matrix copy."""
import argparse,hashlib,json,time
from pathlib import Path
import requests
from bs4 import BeautifulSoup
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--proxy');a=p.parse_args()
assert str(a.root.resolve()).startswith('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/PRAD/B/')
s=a.root/'source';s.mkdir(exist_ok=True);ses=requests.Session()
if a.proxy:ses.proxies={'http':a.proxy,'https':a.proxy}
manifest=[]
for name,id in [('hirz_scrna_annotations.rds','1kJw4UavRdWU0d_fEpuY0TovXwPVZx4t4'),('hirz_scrna_counts.rds','1Fripgmxtww-QgcwmCQOKpf2lOwrPgKyS')]:
    url='https://drive.google.com/uc?export=download&id='+id;dest=s/name
    for attempt in range(3):
        try:
            if not dest.exists():
                x=ses.get(url,stream=True,timeout=(20,90));x.raise_for_status()
                if 'text/html' in x.headers.get('Content-Type',''):
                    form=BeautifulSoup(x.text,'html.parser').find('form',id='download-form')
                    if form is None:raise RuntimeError('No public download form for '+name)
                    params={i['name']:i.get('value','') for i in form.find_all('input') if i.get('name')}
                    x=ses.get(form['action'],params=params,stream=True,timeout=(20,90));x.raise_for_status()
                assert 'text/html' not in x.headers.get('Content-Type','')
                partial=dest.with_suffix('.part')
                with partial.open('wb') as f:
                    for chunk in x.iter_content(1048576):f.write(chunk)
                partial.rename(dest)
            h=hashlib.sha256()
            with dest.open('rb') as f:
                for b in iter(lambda:f.read(1048576),b''):h.update(b)
            manifest.append(dict(file=name,url=url,sha256=h.hexdigest(),bytes=dest.stat().st_size));print(name,dest.stat().st_size,flush=True);break
        except Exception as e:
            print('attempt',attempt+1,type(e).__name__,str(e)[:200],flush=True)
            if attempt==2:raise
    (s/'hirz_scrna_manifest.json').write_text(json.dumps(manifest,indent=2))
