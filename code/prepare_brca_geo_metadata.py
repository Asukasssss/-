"""Fetch GEO annotations on server165, or relay in memory from a network-enabled host.

python prepare_brca_geo_metadata.py --out /absolute/server/run/source [--relay-host server165]
Raw annotations never written on the relay computer. No expression matrix requested.
"""
import argparse,csv,io,json,hashlib,subprocess
from pathlib import Path
import requests

URL='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37751&targ=gsm&form=text&view=quick'

def convert(text):
    samples={};cur=None
    for line in text.splitlines():
        if line.startswith('^SAMPLE = '):
            cur=line.split(' = ',1)[1]
            assert cur not in samples
            samples[cur]={}
        elif line.startswith('!Sample_characteristics_ch1 = ') and cur:
            value=line.split(' = ',1)[1]
            if ': ' in value:
                k,v=value.split(': ',1);samples[cur][k]=v
    assert len(samples)==108
    keys=list(samples);fields=sorted({k for d in samples.values() for k in d})
    stream=io.StringIO();w=csv.writer(stream,delimiter='\t');w.writerow(['!Sample_geo_accession']+keys)
    for field in fields:w.writerow(['!Sample_characteristics_ch1']+[field+': '+samples[k].get(field,'') for k in keys])
    return stream.getvalue().encode()

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',required=True);ap.add_argument('--relay-host');a=ap.parse_args()
    r=requests.get(URL,timeout=(5,25));r.raise_for_status()
    payloads={'GSE37751_quick.txt':r.content,'geo_metadata_lines.txt':convert(r.text)}
    for name,data in payloads.items():
        dest=str(Path(a.out)/name).replace('\\','/')
        if a.relay_host:
            import shlex
            subprocess.run(['ssh',a.relay_host,'test ! -e '+shlex.quote(dest)+' && cat > '+shlex.quote(dest)],input=data,check=True)
        else:
            with open(dest,'xb') as f:f.write(data)
    print(json.dumps({'url':URL,'sha256':{k:hashlib.sha256(v).hexdigest() for k,v in payloads.items()}}))
