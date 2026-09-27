"""Retrieve public GEO metadata and essential Visium files on server165 only."""
from pathlib import Path
import requests, re, json, hashlib, time, concurrent.futures
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/BRCA_GSE210616')
ROOT.mkdir(parents=True,exist_ok=True)
def fetch(url,p):
    if p.exists() and p.stat().st_size>0:return
    for attempt in range(4):
        try:
            with requests.get(url,stream=True,timeout=(30,90)) as r:
                r.raise_for_status()
                with p.with_suffix(p.suffix+'.part').open('wb') as f:
                    for b in r.iter_content(1024*1024):f.write(b)
            p.with_suffix(p.suffix+'.part').replace(p);return
        except Exception:
            if attempt==3:raise
            time.sleep(2+attempt*2)
def one(gsm):
    d=ROOT/gsm;d.mkdir(exist_ok=True)
    p=d/'metadata.soft'
    fetch(f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=quick',p)
    s=p.read_text(); assert '^SAMPLE' in s
    title=re.search(r'!Sample_title = (.+)',s).group(1)
    urls=re.findall(r'!Sample_supplementary_file(?:_\d+)? = (\S+)',s)
    urls=[u.replace('ftp://','https://') for u in urls if '.cloupe' not in u]
    rec={'gsm':gsm,'title':title,'characteristics':re.findall(r'!Sample_characteristics_ch1 = (.+)',s),'files':[]}
    for u in urls:
        p=d/u.rsplit('/',1)[1];fetch(u,p)
        rec['files'].append({'url':u,'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    (d/'manifest.json').write_text(json.dumps(rec,indent=2));print(gsm,title,len(urls),'DONE',flush=True)
    return rec
if __name__=='__main__':
    gsms=re.findall(r'!Series_sample_id = (GSM\d+)',(ROOT/'series_metadata.soft').read_text())
    assert len(gsms)==43
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        records=list(ex.map(one,gsms))
    (ROOT/'manifest.json').write_text(json.dumps(records,indent=2))
    print('ALL_43_DONE',flush=True)
