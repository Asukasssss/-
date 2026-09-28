"""Read public HDF5 metadata through HTTP ranges in RAM; never cache source locally."""
import io,requests,h5py,sys,json
class RemoteH5(io.RawIOBase):
 def __init__(self,url):
  self.url=url;self.s=requests.Session();r=self.s.head(url,timeout=30);r.raise_for_status();self.size=int(r.headers['Content-Length']);self.pos=0;self.cache={};self.transferred=0
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,offset,whence=0):
  self.pos=offset if whence==0 else self.pos+offset if whence==1 else self.size+offset;return self.pos
 def readinto(self,b):
  x=self.read(len(b));b[:len(x)]=x;return len(x)
 def read(self,n=-1):
  if n<0:n=self.size-self.pos
  end=min(self.pos+n,self.size);out=[]
  while self.pos<end:
   block=self.pos//65536;start=block*65536
   if block not in self.cache:
    stop=min(start+65535,self.size-1);r=self.s.get(self.url,headers={'Range':f'bytes={start}-{stop}'},timeout=30)
    r.raise_for_status();assert r.status_code==206 and len(r.content)==stop-start+1
    assert r.headers['Content-Range'].startswith(f'bytes {start}-{stop}/')
    self.cache[block]=r.content;self.transferred+=len(r.content)
    if self.transferred>100*1024*1024:raise RuntimeError('Metadata range budget exceeded')
   b=self.cache[block];k=min(end-self.pos,len(b)-(self.pos-start));out.append(b[self.pos-start:self.pos-start+k]);self.pos+=k
  return b''.join(out)
URL='https://cf.10xgenomics.com/samples/spatial-exp/3.1.2/Visium_HD_FF_Human_Breast_Cancer/Visium_HD_FF_Human_Breast_Cancer_feature_slice.h5'
if __name__=='__main__':
 remote=RemoteH5(URL)
 with h5py.File(remote,'r') as f:
  print('KEYS',list(f.keys()),flush=True)
  if 'features' in f:
   g=f['features'];print('GENE_HITS',[(k,[v.decode() for v in g[k][:] if b'LYPLA1' in v or b'ENSG00000120992' in v]) for k in ['id','name'] if k in g],flush=True)
 print('RANGE_BYTES',remote.transferred,flush=True)
