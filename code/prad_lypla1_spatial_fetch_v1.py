"""Download author-public Hirz data directly to a new server165 run directory.
Usage: python prad_lypla1_spatial_fetch_v1.py RUN_ROOT [--proxy URL]
Use only a new run directory. Never download these matrices to a local/Git clone.
"""
import argparse, json, hashlib
from pathlib import Path
import requests
from bs4 import BeautifulSoup
p=argparse.ArgumentParser();p.add_argument('root');p.add_argument('--proxy');args=p.parse_args()
root=Path(args.root)
assert str(root.resolve()).startswith('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/PRAD/B/')
root.mkdir(exist_ok=False);(root/'.running').write_text('PRAD LYPLA1 spatial v1\n')
src=root/'source';src.mkdir()
s=requests.Session()
if args.proxy:s.proxies={'http':args.proxy,'https':args.proxy}
ids={'hirz_annotations.rds':'1lL721_aXt4-pQIfZcoiUxyvrgnaDtGdO',
     'hirz_coordinates.rds':'10iNhgRAUOVy4RjfMkRnVihdpzJaEk8cC',
     'hirz_counts.rds':'18dqg_SHMNKOtRtzSatnkcjuiACdqTC-L'}
manifest=[]
urls={k:'https://drive.google.com/uc?export=download&id='+v for k,v in ids.items()}
urls['supplementary_data2.xlsx']='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-023-36325-2/MediaObjects/41467_2023_36325_MOESM5_ESM.xlsx'
urls['hirz_article.html']='https://www.nature.com/articles/s41467-023-36325-2'
for name,url in urls.items():
    r=s.get(url,stream=True,timeout=(15,60));r.raise_for_status()
    if name.endswith('.rds') and 'text/html' in r.headers.get('Content-Type',''):
        form=BeautifulSoup(r.text,'html.parser').find('form',id='download-form')
        if form is None:raise RuntimeError('Download returned HTML, not an accessible RDS: '+name)
        params={i['name']:i.get('value','') for i in form.find_all('input') if i.get('name')}
        r=s.get(form['action'],params=params,stream=True,timeout=(15,60));r.raise_for_status()
    if name.endswith('.rds'):assert 'text/html' not in r.headers.get('Content-Type','')
    h=hashlib.sha256()
    with (src/name).open('xb') as f:
        for chunk in r.iter_content(1048576):f.write(chunk);h.update(chunk)
    manifest.append(dict(file=name,url=url,sha256=h.hexdigest(),bytes=(src/name).stat().st_size))
    print(name,(src/name).stat().st_size,flush=True)
(src/'download_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
