"""Resumable, range-validated download of original GSA processed matrices."""
from pathlib import Path
import concurrent.futures,hashlib,json,time,urllib.request
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
DATA=ROOT/'data/candidates/CRA001160_author_v1'
URL='https://download.cncb.ac.cn/gsa/CRA001160/other_files/PDAC.tar.gz'
SIZE=718239348
ETAG='"5f3a14ee-2acf7674"'
CHUNK=4*1024*1024
# Pinned to the actual v1 download, not a claimed author-published checksum.
ARCHIVE_SHA='d4f93c765405d31b9e9f445045de27888dbf1e8b75c425635a078e01e0ea7406'
ANNOTATION_SHA='41a37cfd2a12c06e5583b57d23aab5e97abc6f6a138bc27972badd9bb9477a64'
def file_hash(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(CHUNK),b''):h.update(b)
 return h.hexdigest()
def fetch(i):
 start=i*CHUNK;end=min(SIZE,start+CHUNK)-1
 p=DATA/'chunks'/f'{i:04d}.bin'
 if p.exists() and p.stat().st_size==end-start+1:return i
 for attempt in range(6):
  try:
   req=urllib.request.Request(URL,headers={'Range':f'bytes={start}-{end}','If-Range':ETAG})
   with urllib.request.urlopen(req,timeout=60) as r:
    assert r.status==206 and r.headers['Content-Range']==f'bytes {start}-{end}/{SIZE}'
    assert r.headers['ETag']==ETAG
    b=r.read()
   assert len(b)==end-start+1
   p.write_bytes(b);return i
  except Exception:
   if attempt==5:raise
   time.sleep(2+attempt)
def main():
 DATA.mkdir(exist_ok=True,parents=True);(DATA/'chunks').mkdir(exist_ok=True)
 annotation=DATA/'all_celltype.txt'
 if not annotation.exists():
  with urllib.request.urlopen(URL.rsplit('/',1)[0]+'/all_celltype.txt',timeout=60) as r:content=r.read()
  assert len(content)==2101436 and hashlib.sha256(content).hexdigest()==ANNOTATION_SHA
  annotation.write_bytes(content)
 assert file_hash(annotation)==ANNOTATION_SHA,'Changed annotation input'
 if (DATA/'PDAC.tar.gz').exists() and (DATA/'download_manifest.json').exists():
  assert file_hash(DATA/'PDAC.tar.gz')==ARCHIVE_SHA,'Changed archive input'
  print('PINNED_INPUTS_VERIFIED',flush=True);return
 n=(SIZE+CHUNK-1)//CHUNK
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  for k,_ in enumerate(pool.map(fetch,range(n)),1):
   if k%10==0:print(f'{k}/{n} chunks',flush=True)
 sha=hashlib.sha256()
 with (DATA/'PDAC.tar.gz').open('wb') as out:
  for i in range(n):
   b=(DATA/'chunks'/f'{i:04d}.bin').read_bytes();out.write(b);sha.update(b)
 assert (DATA/'PDAC.tar.gz').stat().st_size==SIZE
 assert sha.hexdigest()==ARCHIVE_SHA,'Archive differs from pinned v1 input'
 (DATA/'download_manifest.json').write_text(json.dumps({'url':URL,'size':SIZE,'etag':ETAG,'sha256':sha.hexdigest(),'range_validation':True},indent=2))
 print('DOWNLOAD_DONE',sha.hexdigest(),flush=True)
if __name__=='__main__':main()
