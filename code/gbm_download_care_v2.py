"""Download public processed CARE counts directly to an exclusive server run."""
from pathlib import Path
import argparse,concurrent.futures,hashlib,json,time
import requests,pandas as pd

def main(out,proxy=None,workers=8):
    assert (out/'.running').exists()
    pinned='https://raw.githubusercontent.com/dravishays/GBM-CARE-WT/383b0a320f0993a4529167cdb591dac63e103674/data/'
    for name,url in [('celltype_meta_data_2025_01_08.RDS',pinned+'celltype_meta_data_2025_01_08.RDS'),('spitzer_supptable1.xlsx',pinned+'spitzer_supptable1.xlsx'),('CARE_filelist.txt','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274546/suppl/filelist.txt')]:
        dest=out/'source'/name
        if not dest.exists():
            r=requests.get(url,timeout=(20,90));r.raise_for_status();dest.write_bytes(r.content)
    assert (out/'source/genes_unique.tsv').exists(),'Copy frozen GBM candidate gene list to source first'
    tab=pd.read_csv(out/'source/CARE_filelist.txt',sep='\t')
    tab=tab[tab.iloc[:,0].eq('File')]
    def fetch(row):
        name=row['Name']; size=int(row['Size']); acc=name.split('_')[0]
        url=f'https://ftp.ncbi.nlm.nih.gov/geo/samples/{acc[:-3]}nnn/{acc}/suppl/{name}'
        dest=out/'source'/name
        for attempt in range(4):
            try:
                if not dest.exists() or dest.stat().st_size!=size:
                    part=dest.with_suffix('.part')
                    start=part.stat().st_size if part.exists() else 0
                    if start>=size:part.unlink();start=0
                    headers={'Range':f'bytes={start}-'} if start else {}
                    route={'https':proxy,'http':proxy} if proxy and attempt%2 else None
                    started=time.monotonic()
                    with requests.get(url,headers=headers,proxies=route,stream=True,timeout=(20,45)) as r:
                        r.raise_for_status()
                        append=start>0 and r.status_code==206
                        if append:assert r.headers.get('Content-Range','').startswith(f'bytes {start}-')
                        with part.open('ab' if append else 'wb') as f:
                            for b in r.iter_content(2**18):
                                f.write(b)
                                if time.monotonic()-started>300:raise TimeoutError('bounded transfer retry with byte-range resume')
                    assert part.stat().st_size==size
                    part.replace(dest)
                return dict(file=name,bytes=size,sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),url=url)
            except Exception:
                if attempt==3:raise
                time.sleep(2+attempt)
    records=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for rec in pool.map(fetch,tab.to_dict('records')):
            records.append(rec); print('downloaded',len(records),'of',len(tab),flush=True)
    pd.DataFrame(records).to_csv(out/'private/CARE_source_file_manifest.tsv',sep='\t',index=False)
    (out/'source/download_complete.json').write_text(json.dumps(dict(files=len(records),bytes=sum(x['bytes'] for x in records))))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--proxy');p.add_argument('--workers',type=int,default=8);a=p.parse_args();main(a.out,a.proxy,a.workers)
