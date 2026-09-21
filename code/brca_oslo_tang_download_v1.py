"""Download processed sources only to isolated server run; record failures."""
from pathlib import Path
import requests,sys,json,hashlib,concurrent.futures,time,subprocess
root=Path(sys.argv[1]);src=root/'source'
files={
'Oslo2_tables.xlsx':'https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs40170-016-0152-x/MediaObjects/40170_2016_152_MOESM2_ESM.xlsx',
'Tang_S2.xlsx':'https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13058-014-0415-9/MediaObjects/13058_2014_415_MOESM2_ESM.xlsx',
'Tang_S1.docx':'https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13058-014-0415-9/MediaObjects/13058_2014_415_MOESM1_ESM.docx',
'GSE58212_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58212/matrix/GSE58212_series_matrix.txt.gz',
'GSE42568_series_matrix.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42568/matrix/GSE42568_series_matrix.txt.gz',
}
def one(item):
 name,url=item;p=src/name
 for attempt in range(3):
  try:
   if not p.exists():
    subprocess.run([sys.executable,'-c',"import requests,sys;from pathlib import Path;x=requests.get(sys.argv[1],timeout=60);x.raise_for_status();assert len(x.content)>1000;Path(sys.argv[2]).write_bytes(x.content)",url,str(p)],timeout=80,check=True,capture_output=True)
   return dict(file=name,url=url,status='DONE',bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
  except Exception as e:
   if attempt==2:return dict(file=name,url=url,status='ACCESS_BLOCKED',reason=str(e)[:200])
   time.sleep(1)
rows=list(concurrent.futures.ThreadPoolExecutor(max_workers=3).map(one,files.items()))
(root/'public/download_manifest.json').write_text(json.dumps(rows,indent=2)+'\n')
for r in rows:print(r,flush=True)
