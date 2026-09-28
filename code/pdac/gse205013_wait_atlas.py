"""Server job: validate atlas download, then inspect selected-study labels."""
from pathlib import Path
import time,json,hashlib,gzip,subprocess,shutil
R=Path(__file__).resolve().parent
D=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/PDAC_Loveless_14199536')
f=D/'scAtlas.rds.gz'
def state(status,reason):
 (R/'public/atlas_job_status.json').write_text(json.dumps(dict(status=status,reason=reason,updated_unix=time.time()),indent=2));print(status,reason,flush=True)
def main():
 state('PARTIAL','Waiting for atlas download; malignancy comparisons NOT_RUN')
 deadline=time.time()+48*3600
 while (D/'scAtlas.rds.gz.aria2').exists() or not f.exists():
  running=False
  for proc in Path('/proc').glob('[0-9]*/cmdline'):
   try:
    cmd=proc.read_bytes().replace(b'\x00',b' ')
    if b'aria2c' in cmd and b'-o scAtlas.rds.gz ' in cmd:running=True
   except (FileNotFoundError,PermissionError,ProcessLookupError):pass
  if not running:raise RuntimeError('Atlas downloader stopped before completion; retain partial download for resume')
  if time.time()>deadline:raise RuntimeError('48h download wait limit reached; extraction not run')
  time.sleep(30)
 assert f.stat().st_size==33413868884
 state('PARTIAL','Download complete; verifying original archive MD5/SHA256')
 md5=hashlib.md5();sha=hashlib.sha256()
 with f.open('rb') as inp:
  for block in iter(lambda:inp.read(8*1024*1024),b''):md5.update(block);sha.update(block)
 assert md5.hexdigest()=='705078352e1feb260a64cec67c64ade0','Zenodo MD5 mismatch'
 (R/'public/atlas_source_verified.json').write_text(json.dumps(dict(url='https://zenodo.org/records/14199536/files/scAtlas.rds.gz',bytes=f.stat().st_size,md5=md5.hexdigest(),sha256=sha.hexdigest()),indent=2))
 with gzip.open(f,'rb') as inp:magic=inp.read(2)
 source=f
 if magic==b'\x1f\x8b':
  state('PARTIAL','Verified double gzip; unpacking outer layer on server')
  source=D/'scAtlas.inner.rds';part=D/'scAtlas.inner.rds.part'
  with gzip.open(f,'rb') as inp,part.open('wb') as out:shutil.copyfileobj(inp,out,8*1024*1024)
  part.rename(source)
 mem={line.split(':')[0]:int(line.split()[1]) for line in Path('/proc/meminfo').read_text().splitlines()}
 if mem['MemAvailable']<260*1024*1024:raise RuntimeError('Less than260GiB available; defer memory-intensive atlas loading')
 state('PARTIAL','Reading verified atlas; export selected GSE205013 metadata and LYPLA1 on server only')
 with (R/'atlas_extract.log').open('w') as log:
  subprocess.run(['Rscript',str(R/'gse205013_atlas_schema.R'),str(source),str(R)],stdout=log,stderr=subprocess.STDOUT,check=True)
 state('NEEDS_REVIEW','Atlas selected-study extraction complete; review label identity and match GEO counts before contrasts')
if __name__=='__main__':
 try:main()
 except Exception as e:state('ACCESS_BLOCKED',str(e));raise
