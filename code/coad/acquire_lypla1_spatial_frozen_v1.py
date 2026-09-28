"""Acquire all 12 GSE283052 frozen CRC sections on server165.

Run before analyze_lypla1_spatial_v1.py. A memory-only HTTP-to-SFTP relay
may prepopulate complete files when server direct egress is slow.
"""
import concurrent.futures,gzip,json,urllib.request
from pathlib import Path
DATA=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/coad_lypla1_spatial_20260928')
assert DATA.is_dir()
softurl='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/GSE283052/soft/GSE283052_family.soft.gz'

def fetch(url):
    path=DATA/url.rsplit('/',1)[1]
    if path.exists():return str(path)
    with urllib.request.urlopen(url,timeout=90) as response,path.with_suffix(path.suffix+'.part').open('wb') as stream:
        expected=int(response.headers.get('Content-Length',0));n=0
        while True:
            block=response.read(1024*1024)
            if not block:break
            stream.write(block);n+=len(block)
        assert not expected or n==expected
    path.with_suffix(path.suffix+'.part').replace(path)
    print('READY',path.name,n,flush=True)
    return str(path)

soft=gzip.open(fetch(softurl),'rt').read()
urls=[x.split(' = ',1)[1].replace('ftp://','https://',1) for x in soft.splitlines()
    if x.startswith(('!Sample_supplementary_file','!Series_supplementary_file'))]
suffixes=['matrix.mtx.gz','features.tsv.gz','barcodes.tsv.gz','scalefactors_json.json.gz','tissue_lowres_image.png.gz','tissue_positions_list.csv.gz']
urls=sorted({u for u in urls if any(u.endswith(s) for s in suffixes)})
assert len(urls)==50 and sum(u.endswith('matrix.mtx.gz') for u in urls)==12
(DATA/'frozen_acquisition_lock.json').write_text(json.dumps(dict(
    cohort='GSE283052',scope='all 12 CRC sections before expression results',urls=urls),indent=2)+'\n')
with concurrent.futures.ThreadPoolExecutor(3) as pool:
    list(pool.map(fetch,urls))
print('FROZEN_ACQUIRED',flush=True)
