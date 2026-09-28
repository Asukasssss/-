"""Record an actually verified delivery commit, then refresh package integrity."""
from pathlib import Path
import json,hashlib,zipfile,subprocess,sys
import pandas as pd
from complete_report_v5 import ROOT,OUT,PUBLIC,delivery_files,tsv,read,normalize_public
def main():
 head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
 remote=subprocess.check_output(['git','-c','http.version=HTTP/1.1','ls-remote','origin','refs/heads/analysis/pdac-initial'],cwd=ROOT,text=True).split()[0]
 assert head==remote
 evidence=dict(status='DONE',verified_delivery_commit=head,branch='analysis/pdac-initial',method='git ls-remote exact SHA equality',PR=2,draft=True,merged=False,note='This records the published report commit; the later record-only commit is not self-referenced.')
 for folder in [OUT,PUBLIC]:
  (folder/'publication_verified.json').write_text(json.dumps(evidence,indent=2)+'\n',encoding='utf8')
  d=read(folder/'acceptance.tsv');d.loc[d.check_id.eq('14'),['status','evidence_path','reason']]=['DONE','publication_verified.json','Actual public report commit verified remotely; PR remains draft and unmerged']
  tsv(folder/'acceptance.tsv',d)
 for filename in ['prepare_report_v5.py','configure_report_v5.py','render_report_v5.py','build_report_v5_workbook.mjs','complete_report_v5.py','verify_report_v5.py','refine_report_v5_layout.py','register_report_v5.py','record_report_v5_publication.py']:
  for folder in [OUT,PUBLIC]:(folder/'repro'/filename).write_bytes((ROOT/'code/pdac'/filename).read_bytes())
 files=sorted(delivery_files());checks=pd.DataFrame([dict(path=p.relative_to(OUT).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]);tsv(OUT/'checksums.tsv',checks)
 zpath=OUT.parent/'PDAC_统一报告完整包_v5.zip'
 with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
  for p in files+[OUT/'checksums.tsv']:z.write(p,'PDAC_v5/'+p.relative_to(OUT).as_posix())
 with zipfile.ZipFile(zpath) as z:
  assert z.testzip() is None
  for r in checks.itertuples():assert hashlib.sha256(z.read('PDAC_v5/'+r.path)).hexdigest()==r.sha256
 info=dict(status='PASS',files=len(files)+1,zip_sha256=hashlib.sha256(zpath.read_bytes()).hexdigest(),zip_bytes=zpath.stat().st_size)
 for folder in [OUT,PUBLIC]:(folder/'package_validation.json').write_text(json.dumps(info,indent=2)+'\n',encoding='utf8')
 pubs=sorted(p for p in PUBLIC.rglob('*') if p.is_file() and p.name!='checksums.tsv');tsv(PUBLIC/'checksums.tsv',pd.DataFrame([dict(path=p.relative_to(PUBLIC).as_posix(),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in pubs]))
 normalize_public()
 print(json.dumps(dict(evidence,package=info)))
if __name__=='__main__':main()
