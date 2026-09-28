"""Acquire author-selected one-section-per-case Space Ranger archives on server165."""
import concurrent.futures, hashlib, json, time, urllib.request
from pathlib import Path
DATA=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/coad_lypla1_pathology_20260928')
SAMPLES=[
 ('S1_Cec_Rep1','SN123_A551763_Rep1','b26940f8bf3b3e9855b0a116825e638c'),
 ('S2_Col_R_Rep1','SN123_A595688_Rep1','b5aa2ee18977b0be67a72ff2b966ab26'),
 ('S3_Col_R_Rep1','SN048_A416371_Rep1','8aec27998074f672fcdaba6faa609da0'),
 ('S4_Col_Sig_Rep1','SN84_A120838_Rep1','4a185c5e1bf88995ee0c03dbb4e8b63c'),
 ('S5_Rec_Rep1','SN048_A121573_Rep1','608a39f21da059024121407967a76b8a'),
 ('S6_Rec_Rep2','SN124_A938797_Rep2','ee145e0f97146586baa8ef95edae05b0'),
 ('S7_RecSig_Rep1','SN123_A798015_Rep1','0e589a2c96546fff5ca21cb18d597791')]
def acquire(row):
 alias,name,md5=row; p=DATA/(name+'.zip')
 if p.exists() and hashlib.md5(p.read_bytes()).hexdigest()==md5:
  print('REUSE',alias,flush=True);return
 url='https://zenodo.org/api/records/7760264/files/'+p.name+'/content'
 for attempt in range(3):
  try:
   h=hashlib.md5()
   with urllib.request.urlopen(url,timeout=60) as r, p.with_suffix('.part').open('wb') as f:
    while True:
     b=r.read(1024*1024)
     if not b:break
     f.write(b);h.update(b)
   assert h.hexdigest()==md5,(name,h.hexdigest())
   p.with_suffix('.part').replace(p)
   print('DONE',alias,p.stat().st_size,flush=True);return
  except Exception as e:
   print('RETRY',alias,str(e),flush=True)
   if attempt==2:raise
   time.sleep(3)
if __name__=='__main__':
 DATA.mkdir(exist_ok=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
  list(pool.map(acquire,SAMPLES))
 print('ACQUISITION_DONE',flush=True)
