"""Stream public supplements directly into SFTP; no local source files."""
import sys,getpass,urllib.request,json,hashlib
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'runtime/ssh_dependencies'))
import paramiko
c=paramiko.SSHClient();c.load_host_keys(str(Path.home()/'.ssh/known_hosts'))
p=getpass.getpass('SSH password: ');c.connect('172.22.148.165',username='xuzx',password=p,look_for_keys=False,allow_agent=False);del p
base='/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/'
run=base+'results/collaborative/COAD/B/20260921T050501Z_integral_acquire_v1/'
dest=base+'data/candidates/coad_integral_20260921/'
records=[]
with c.open_sftp()as s:
 for aid in [28156574,28156568]:
  meta=json.load(urllib.request.urlopen('https://api.figshare.com/v2/articles/'+str(aid),timeout=30))
  assert meta['resource_doi']=='10.1021/acs.analchem.4c04421'
  with s.open(dest+str(aid)+'.json','w')as f:f.write(json.dumps(meta))
  for x in meta['files']:
   path=dest+x['name'];md5=hashlib.md5();sha=hashlib.sha256();size=0
   with urllib.request.urlopen(x['download_url'],timeout=120)as h,s.open(path,'wx')as f:
    while True:
     b=h.read(1024*1024)
     if not b:break
     f.write(b);md5.update(b);sha.update(b);size+=len(b)
   assert size==x['size'] and md5.hexdigest()==x['computed_md5']
   records.append(dict(article=aid,doi=meta['doi'],url=x['download_url'],path=path,bytes=size,md5=md5.hexdigest(),sha256=sha.hexdigest(),transfer='HTTPS memory chunks directly to server165 SFTP; no local source file'))
   print('VERIFIED',x['name'],size,flush=True)
 with s.open(run+'source_files.json','w')as f:f.write(json.dumps(records,indent=2))
c.close()
