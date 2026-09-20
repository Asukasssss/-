"""Download public cell-source inputs only to server165; never export patient rows."""
import argparse,concurrent.futures,hashlib,json,urllib.request
from pathlib import Path
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
def main():
 p=argparse.ArgumentParser();p.add_argument('--data-dir',type=Path,required=True);p.add_argument('--mode',choices=['metadata','matrices'],required=True);a=p.parse_args();d=a.data_dir.resolve();assert d.parent==ROOT/'data/candidates';d.mkdir(exist_ok=True)
 geo='https://ftp.ncbi.nlm.nih.gov/geo/series/'
 urls={
 'lee_annotation.txt.gz':geo+'GSE132nnn/GSE132465/suppl/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz',
 'pelka_cluster.csv.gz':geo+'GSE178nnn/GSE178341/suppl/GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz',
 'pelka_metadata.csv.gz':geo+'GSE178nnn/GSE178341/suppl/GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz'}
 for acc in ['GSE132465','GSE132257','GSE144735','GSE178341','GSE89076']:urls[acc+'_family.soft.gz']=geo+acc[:-3]+'nnn/'+acc+'/soft/'+acc+'_family.soft.gz'
 if a.mode=='matrices':urls={'lee_raw_UMI.txt.gz':geo+'GSE132nnn/GSE132465/suppl/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz','pelka_raw.h5':geo+'GSE178nnn/GSE178341/suppl/GSE178341_crc10x_full_c295v4_submit.h5'}
 def fetch(item):
  name,u=item;dest=d/name
  if dest.exists():raise FileExistsError(str(dest))
  try:
   with urllib.request.urlopen(u,timeout=60) as r,dest.with_suffix(dest.suffix+'.part').open('xb') as f:
    while True:
     b=r.read(1024*1024)
     if not b:break
     f.write(b)
   dest.with_suffix(dest.suffix+'.part').rename(dest)
   h=hashlib.sha256()
   with dest.open('rb') as f:
    for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
   result={'name':name,'url':u,'status':'DONE','bytes':dest.stat().st_size,'sha256':h.hexdigest()}
  except Exception as e:result={'name':name,'url':u,'status':'ACCESS_BLOCKED','error':str(e)}
  print(json.dumps(result),flush=True);return result
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:results=list(ex.map(fetch,urls.items()))
 with (d/('download_'+a.mode+'.json')).open('x') as f:json.dump(results,f,indent=2)
if __name__=='__main__':main()
