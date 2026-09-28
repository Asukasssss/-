"""HTTP Range select required ZIP members in memory, CRC verify, write only server165.

No source data is written locally. Whole source ZIP MD5 is NOT verified in this mode.
"""
import concurrent.futures, getpass, hashlib, io, json, struct, sys, time, urllib.request, zipfile, zlib
from pathlib import Path
from acquire_lypla1_pathology_v1 import SAMPLES, DATA
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'runtime/ssh_dependencies'))
import paramiko

def fetch(u,start,end):
 for k in range(4):
  try:
   req=urllib.request.Request(u,headers={'Range':f'bytes={start}-{end}'})
   with urllib.request.urlopen(req,timeout=40) as r:
    assert r.status==206
    assert r.headers['Content-Range'].startswith(f'bytes {start}-')
    b=r.read();assert len(b)==end-start+1
   return b
  except Exception as err:
   print("RANGE_RETRY",start,end,repr(err),flush=True)
   if k==3:raise
   time.sleep(2+k)

def main():
 c=paramiko.SSHClient();c.load_host_keys(str(Path.home()/'.ssh/known_hosts'))
 c.connect('172.22.148.165',username='xuzx',password=getpass.getpass('SSH password (not saved): '),look_for_keys=False,allow_agent=False,timeout=20)
 dest=DATA.as_posix()+'/'
 print('SSH_READY',flush=True)
 for alias,name,md5 in SAMPLES:
  print('START',alias,flush=True)
  with c.open_sftp() as s:
   try:s.stat(dest+name+'.selected.zip');print('EXISTS',alias,flush=True);continue
   except FileNotFoundError:pass
  u='https://zenodo.org/api/records/7760264/files/'+name+'.zip/content'
  with urllib.request.urlopen(urllib.request.Request(u,headers={'Range':'bytes=-65536'}),timeout=40) as r:
   tail=r.read();size=int(r.headers['Content-Range'].split('/')[-1])
  i=tail.rfind(b'PK\x05\x06');e=list(struct.unpack('<4s4H2LH',tail[i:i+22]));cdsize,offset=e[5:7]
  assert e[-1]==0
  cd=fetch(u,offset,offset+cdsize-1);e[6]=0
  z=zipfile.ZipFile(io.BytesIO(cd+struct.pack('<4s4H2LH',*e)))
  suffixes=['filtered_feature_bc_matrix.h5','tissue_positions_list.csv','scalefactors_json.json','tissue_lowres_image.png']
  selected=[]
  for suffix in suffixes:
   match=[v for v in z.infolist() if v.filename.endswith(suffix) and '__MACOSX' not in v.filename];assert len(match)==1
   selected+=match
  buf=io.BytesIO();proof=[]
  with zipfile.ZipFile(buf,'w',compression=zipfile.ZIP_DEFLATED) as out:
   for v in selected:
    header=fetch(u,v.header_offset,v.header_offset+29);h=struct.unpack('<4s5H3L2H',header);assert h[0]==b'PK\x03\x04'
    start=v.header_offset+30+h[-2]+h[-1];end=start+v.compress_size
    print('FETCH',v.filename,start,end,flush=True)
    spans=[(p,min(p+262144,end)-1) for p in range(start,end,262144)]
    with concurrent.futures.ThreadPoolExecutor(4) as pool:
     blocks=list(pool.map(lambda ab:fetch(u,*ab),spans))
    compressed=b''.join(blocks)
    data=zlib.decompress(compressed,-15) if v.compress_type==8 else compressed
    assert len(data)==v.file_size and zlib.crc32(data)&0xffffffff==v.CRC
    out.writestr(v.filename,data)
    proof.append(dict(member=v.filename,bytes=len(data),crc32=v.CRC,sha256=hashlib.sha256(data).hexdigest()))
    print('MEMBER',alias,v.filename,len(data),flush=True)
  b=buf.getvalue();manifest=dict(source_url=u,source_zip_bytes=size,official_source_md5=md5,source_whole_md5_verified=False,member_crc_verified=True,selected_zip_sha256=hashlib.sha256(b).hexdigest(),members=proof)
  with c.open_sftp() as s:
   s.putfo(io.BytesIO(b),dest+name+'.selected.zip')
   s.putfo(io.BytesIO((json.dumps(manifest,indent=2)+'\n').encode()),dest+name+'.selected_manifest.json')
  print('READY',alias,len(b),flush=True)
 c.close();print('ALL_READY',flush=True)
if __name__=='__main__':main()
