import concurrent.futures as cf,subprocess,re,gzip,hashlib,json,time
from pathlib import Path
s=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/ccrcc_lypla1_gse269826_20260928');u='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269826/suppl/GSE269826_Freshbiopsies_epithelial_reclustered.rds.gz'
chunks=s/'verified_ranges';chunks.mkdir(exist_ok=True)
x=subprocess.run(['curl','-4','-sS','-L','--range','0-0','--max-time','30','-D','/dev/stderr',u],capture_output=True);total=int(re.search(rb'content-range:\s*bytes 0-0/(\d+)',x.stderr,re.I)[1]);block=8388608
state=dict(url=u,total_bytes=total,range_bytes=block,workers=6,status='DOWNLOADING',completed_bytes=0)
def save(): (s/'range_status.json').write_text(json.dumps(state,indent=2))
def one(a):
 b=min(a+block,total)-1;p=chunks/('%d-%d.bin'%(a,b))
 if p.exists() and p.stat().st_size==b-a+1:return p
 for k in range(4):
  x=subprocess.run(['curl','-4','-sS','-L','--fail','--range','%d-%d'%(a,b),'--connect-timeout','15','--max-time','90','-D','/dev/stderr',u],capture_output=True)
  if x.returncode==0 and re.search(('content-range:\\s*bytes %d-%d/%d'%(a,b,total)).encode(),x.stderr,re.I) and len(x.stdout)==b-a+1:
   p.write_bytes(x.stdout);return p
 raise RuntimeError('Failed checked range %d'%a)
save()
with cf.ThreadPoolExecutor(max_workers=6) as ex:
 for f in cf.as_completed([ex.submit(one,a) for a in range(0,total,block)]):
  p=f.result();state['completed_bytes']+=p.stat().st_size;save()
state['status']='VERIFYING';save();dest=s/'epithelial.rds.gz';h=hashlib.sha256()
with dest.with_suffix('.assembling').open('wb') as f:
 for a in range(0,total,block):
  z=(chunks/('%d-%d.bin'%(a,min(a+block,total)-1))).read_bytes();f.write(z);h.update(z)
dest.with_suffix('.assembling').replace(dest)
with gzip.open(dest,'rb') as f, (s/'epithelial_inner.rds').open('wb') as g:
 for z in iter(lambda:f.read(1048576),b''):g.write(z)
state.update(status='COMPLETE',sha256=h.hexdigest(),gzip_crc_verified=True);save();print(json.dumps(state),flush=True)


