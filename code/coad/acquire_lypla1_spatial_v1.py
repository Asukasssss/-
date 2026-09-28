"""Acquire selected original GEO Visium HD outputs on server165 only."""
import argparse,concurrent.futures,gzip,hashlib,json,re,time,urllib.request
from pathlib import Path
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
DATA=ROOT/'data/candidates/coad_lypla1_spatial_20260928'
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--commit',required=True);a=ap.parse_args()
assert a.out.parent==ROOT/'results/collaborative/COAD/B' and (a.out/'.running').is_dir()
DATA.mkdir(exist_ok=True);pub=a.out/'public';pub.mkdir(exist_ok=True)
samples={'GSM8594567':'P1CRC','GSM8594568':'P2CRC','GSM8594569':'P5CRC'}
suffixes=['filtered_feature_bc_matrix.h5','Metadata.parquet.gz','tissue_positions.parquet.gz','scalefactors_json.json.gz','tissue_lowres_image.png.gz','probe_set.csv.gz']
jobs=[]
for gsm,sample in samples.items():
    base='https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8594nnn/'+gsm+'/suppl/'
    listing=urllib.request.urlopen(base,timeout=40).read().decode()
    (DATA/(gsm+'_listing.html')).write_text(listing)
    for suffix in suffixes:
        name=gsm+'_'+sample+'_'+suffix;assert 'href="'+name+'"' in listing
        jobs.append((sample,base+name,DATA/name))
jobs.append(('study','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE280nnn/GSE280315/soft/GSE280315_family.soft.gz',DATA/'GSE280315_family.soft.gz'))
def download(job):
    sample,url,path=job
    if not path.exists():
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url,timeout=60) as r,path.with_suffix(path.suffix+'.part').open('wb') as f:
                    length=int(r.headers.get('Content-Length',0));n=0
                    while True:
                        b=r.read(1024*1024)
                        if not b:break
                        f.write(b);n+=len(b)
                    assert not length or n==length
                path.with_suffix(path.suffix+'.part').replace(path);break
            except Exception:
                if attempt==2:raise
                time.sleep(2)
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    print('READY',path.name,path.stat().st_size,flush=True)
    return dict(sample=sample,source_id=path.name,url=url,server_path=str(path),bytes=path.stat().st_size,sha256=h.hexdigest())
with concurrent.futures.ThreadPoolExecutor(3) as pool:manifest=list(pool.map(download,jobs))
(pub/'acquisition_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
(pub/'acquisition_spec.json').write_text(json.dumps(dict(code_lock_commit=a.commit,cohort='GSE280315',samples=samples,selection='all three available tumor Visium HD samples;before examining LYPLA1 values',resolution='8um bins to verify from actual barcodes',question='LYPLA1 spatial expression and author region context',no_new_tests=True,raw_data_stay_server=True),indent=2)+'\n')
(a.out/'ACQUIRED').write_text('DONE\n')
print('ACQUIRED',flush=True)
