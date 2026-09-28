"""Extract selected public ZIP members through RAM directly to server165."""
import io,zipfile,requests,hashlib,json,datetime,sys,concurrent.futures,struct,zlib,threading,time
from brca_gse210616_stream_to_server import remote,transfer
ROOT='/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716'
D=ROOT+'/data/candidates/BRCA_CTA2025'
RUN='20260928T060000Z_cta_lypla1_v1'
OUT=ROOT+'/results/collaborative/BRCA/A/'+RUN
class RangeFile(io.RawIOBase):
 def __init__(self,url):
  self.url=url;self.s=requests.Session();r=self.s.head(url,timeout=45);r.raise_for_status();self.size=int(r.headers['Content-Length']);self.pos=0;self.cache={}
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,o,w=0):
  self.pos=o if w==0 else self.pos+o if w==1 else self.size+o;return self.pos
 def read(self,n=-1):
  n=min(n if n>=0 else self.size-self.pos,self.size-self.pos)
  if not n:return b''
  end=self.pos+n;out=[];bs=65536
  while self.pos<end:
   start=self.pos//bs*bs;stop=min(start+bs-1,self.size-1)
   if start not in self.cache:
    with self.s.get(self.url,headers={'Range':f'bytes={start}-{stop}'},timeout=30,stream=True) as r:
     r.raise_for_status();assert r.status_code==206 and r.headers['Content-Range'].startswith(f'bytes {start}-{stop}/')
     b=r.content;assert len(b)==stop-start+1
    if len(self.cache)>512:self.cache.clear()
    self.cache[start]=b
   b=self.cache[start];k=min(end-self.pos,len(b)-(self.pos-start));out.append(b[self.pos-start:self.pos-start+k]);self.pos+=k
  return b''.join(out)
def main():
 if '--resume-own-run' not in sys.argv:remote(f"mkdir -p '{D}' '{OUT}'; mkdir '{OUT}/.running'")
 records=[];jobs=[]
 if '--resume-own-run' in sys.argv:
  try:records=json.loads(remote("cat '"+D+"/manifest.json'"))
  except Exception:pass
 completed={r['file'] for r in records}
 for archive in ['spaceranger_output.zip','Images.zip']:
  url='https://zenodo.org/api/records/15211538/files/'+archive+'/content'
  z=zipfile.ZipFile(RangeFile(url));infos=z.infolist()
  inventory=[{'name':i.filename,'size':i.file_size,'compressed':i.compress_size} for i in infos]
  transfer([json.dumps(inventory).encode()],D+'/'+archive+'.inventory.json')
  for i in infos:
   n=i.filename
   take=('README' in n or n.endswith('/filtered_feature_bc_matrix.h5') or ('/spatial/' in n and n.endswith(('tissue_lowres_image.png','tissue_positions_list.csv','scalefactors_json.json')))) if archive.startswith('spaceranger') else ('/Manual_annotation/' in n)
   if i.is_dir() or not take or n in completed:continue
   jobs.append((url,i))
 lock=threading.Lock()
 def one(job):
  url,i=job;n=i.filename;path=D+'/'+n
  # Header + compressed member; bounded request validated before body consumption.
  def get(a,b):
   for attempt in range(3):
    try:
     with requests.get(url,headers={'Range':f'bytes={a}-{b}'},stream=True,timeout=60) as r:
      r.raise_for_status();assert r.status_code==206 and r.headers['Content-Range'].startswith(f'bytes {a}-{b}/')
      v=r.content;assert len(v)==b-a+1;return v
    except Exception:
     if attempt==2:raise
  head=get(i.header_offset,i.header_offset+29);assert head[:4]==b'PK\x03\x04'
  nl,el=struct.unpack_from('<HH',head,26);start=i.header_offset+30+nl+el
  compressed=get(start,start+i.compress_size-1) if i.compress_size else b''
  content=zlib.decompress(compressed,-15) if i.compress_type==8 else compressed
  assert len(content)==i.file_size and zlib.crc32(content)==i.CRC
  with lock:
   for attempt in range(4):
    try:
     remote("mkdir -p '"+path.rsplit('/',1)[0]+"'");transfer([content],path);break
    except Exception:
     if attempt==3:raise
     time.sleep(1)
   records.append({'file':n,'url':url,'member':n,'sha256':hashlib.sha256(content).hexdigest(),'bytes':i.file_size})
   for attempt in range(4):
    try:transfer([json.dumps(records,indent=2).encode()],D+'/manifest.json');break
    except Exception:
     if attempt==3:raise
     time.sleep(1)
  print(n,'DONE',flush=True)
 print('REMAINING',len(jobs),flush=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=12) as ex:
  futures={ex.submit(one,j):j[1].filename for j in sorted(jobs,key=lambda j:j[1].file_size)}
  for f in concurrent.futures.as_completed(futures):
   try:f.result()
   except Exception as e:print('FAILED',futures[f],repr(e),flush=True)
 print('FETCH_DONE',len(records),flush=True)
if __name__=='__main__':main()
