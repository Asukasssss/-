"""Public Erickson2022 v5 annotated Visium sections; download on server165 only.
Usage: python prad_lypla1_erickson_fetch_v1.py RUN_ROOT [--proxy URL]
Raw counts and spot annotations must never be copied to the Git workspace.
"""
import argparse,json,hashlib,time,concurrent.futures
from pathlib import Path
import requests
p=argparse.ArgumentParser();p.add_argument('root');p.add_argument('--proxy');args=p.parse_args()
root=Path(args.root)
assert str(root.resolve()).startswith('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/PRAD/B/')
root.mkdir(exist_ok=True);src=root/'source';src.mkdir(exist_ok=True)
lock=root/'.running'
if lock.exists():assert lock.read_text().strip()==root.name
else:lock.write_text(root.name)
proxy={'https':args.proxy} if args.proxy else None
base='https://data.mendeley.com/public-api/datasets/svw96g68dv'
def req(url,**kwargs):
    for k in range(4):
        try:
            r=requests.get(url,proxies=proxy,timeout=(15,45),**kwargs);r.raise_for_status();return r
        except requests.RequestException:
            if k==3:raise
            time.sleep(1+k)
folders=req(base+'/folders/5').json();(src/'folders.json').write_text(json.dumps(folders,indent=2))
by={x['id']:x for x in folders}
def path(x):return (path(by[x['parent_id']])+'/' if x.get('parent_id') else '')+x['name']
selected=[x for x in folders if path(x).startswith('Count_matrices/Patient 1/Visium_with_annotation/') or path(x)=='10x_spaceranger_tissue_highres_images/Patient 1']
files=[]
sections={x['name'] for x in selected if path(x).startswith('Count_matrices/')}
assert len(sections)==7
for x in selected:
    rows=req(base+'/files',params={'folder_id':x['id'],'version':5},headers={'Accept':'application/vnd.mendeley-public-dataset.1+json'}).json()
    for a in rows:
        sec=a['filename'].split('_tissue_hires')[0] if 'highres_images' in path(x) else x['name']
        if sec in sections:a['section']=sec;a['folder_path']=path(x);files.append(a)
(src/'selected_inventory.json').write_text(json.dumps(files,indent=2))
def download(a):
    dest=src/'Patient1'/a['section']/a['filename'];dest.parent.mkdir(parents=True,exist_ok=True)
    cd=a['content_details'];expected=cd['sha256_hash']
    if not dest.exists() or hashlib.sha256(dest.read_bytes()).hexdigest()!=expected:
        for attempt in range(4):
            try:
                r=req(cd['download_url'],stream=True)
                with dest.open('wb') as f:
                    for b in r.iter_content(1048576):f.write(b)
                assert hashlib.sha256(dest.read_bytes()).hexdigest()==expected
                break
            except (requests.RequestException,AssertionError):
                if attempt==3:raise
    digest=hashlib.sha256(dest.read_bytes()).hexdigest();assert digest==expected,(dest,digest,expected)
    print(a['section'],a['filename'],'verified',flush=True)
    return dict(file=str(dest.relative_to(root)),url=cd['download_url'],sha256=digest,bytes=dest.stat().st_size,source_version=5)
with concurrent.futures.ThreadPoolExecutor(3) as ex:manifest=list(ex.map(download,files))
(src/'download_manifest.json').write_text(json.dumps(manifest,indent=2))
print('All selected files verified:',len(manifest))
