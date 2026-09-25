"""Download public processed CARE counts directly to an exclusive server run."""
from pathlib import Path
import argparse,concurrent.futures,hashlib,json,time
import requests,pandas as pd

def main(out):
    assert (out/'.running').exists()
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
                    with requests.get(url,stream=True,timeout=(20,90)) as r:
                        r.raise_for_status()
                        with part.open('wb') as f:
                            for b in r.iter_content(2**20):f.write(b)
                    assert part.stat().st_size==size
                    part.replace(dest)
                return dict(file=name,bytes=size,sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),url=url)
            except Exception:
                if attempt==3:raise
                time.sleep(2+attempt)
    records=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        for rec in pool.map(fetch,tab.to_dict('records')):
            records.append(rec); print('downloaded',len(records),'of',len(tab),flush=True)
    pd.DataFrame(records).to_csv(out/'private/CARE_source_file_manifest.tsv',sep='\t',index=False)
    (out/'source/download_complete.json').write_text(json.dumps(dict(files=len(records),bytes=sum(x['bytes'] for x in records))))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);main(p.parse_args().out)
