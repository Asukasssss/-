"""Server-only extraction of LYPLA1 from GEO counts; no malignancy inferred."""
from pathlib import Path
import tarfile,gzip,io,json,hashlib
import numpy as np,pandas as pd
R=Path(__file__).resolve().parent;P=R/'public';V=R/'private';P.mkdir(exist_ok=True);V.mkdir(exist_ok=True)
SOURCE=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/GSE205013/GSE205013_RAW.tar')
def main():
 rows=[];total=0
 with tarfile.open(SOURCE) as t:
  members={m.name:m for m in t.getmembers() if m.isfile()};mats=sorted(n for n in members if n.endswith('_matrix.mtx.gz'));assert len(mats)==27
  assert len(members)==81
  for name in mats:
   prefix=name[:-len('_matrix.mtx.gz')];sample=prefix.split('_')[-1]
   features=pd.read_csv(gzip.GzipFile(fileobj=t.extractfile(members[prefix+'_features.tsv.gz'])),sep='\t',header=None)
   barcodes=pd.read_csv(gzip.GzipFile(fileobj=t.extractfile(members[prefix+'_barcodes.tsv.gz'])),sep='\t',header=None)[0]
   assert barcodes.notna().all() and barcodes.is_unique
   ix=np.flatnonzero(features[1].eq('LYPLA1').to_numpy());assert len(ix)==1;gene=int(ix[0])+1
   with gzip.GzipFile(fileobj=t.extractfile(members[name])) as f:
    line=f.readline()
    while line.startswith(b'%'):line=f.readline()
    ng,nc,nnz=map(int,line.split());assert ng==len(features) and nc==len(barcodes)
    lib=np.zeros(nc,np.int64);ly=np.zeros(nc,np.int64);seen=0
    while True:
     chunk=f.read(4*1024*1024)
     if not chunk:break
     chunk+=f.readline();z=np.fromstring(chunk.decode('ascii'),sep=' ',dtype=np.int64).reshape(-1,3)
     assert (z[:,0]>=1).all() and (z[:,0]<=ng).all() and (z[:,1]>=1).all() and (z[:,1]<=nc).all() and (z[:,2]>0).all()
     np.add.at(lib,z[:,1]-1,z[:,2]);a=z[z[:,0]==gene];np.add.at(ly,a[:,1]-1,a[:,2]);seen+=len(z)
    assert seen==nnz and (lib>0).all()
   out=pd.DataFrame(dict(sample=sample,barcode=barcodes,LYPLA1_count=ly,total_UMI=lib,LYPLA1_log1p=np.log1p(10000*ly/lib)))
   out.to_csv(V/(sample+'_LYPLA1.tsv.gz'),sep='\t',index=False)
   rows.append(dict(sample=sample,n_barcodes=nc,n_genes=ng,n_entries=nnz,LYPLA1_unique=True));total+=nc;print('EXTRACTED',sample,nc,flush=True)
 pd.DataFrame(rows).to_csv(V/'sample_input_inventory.tsv',sep='\t',index=False)
 h=hashlib.sha256()
 with SOURCE.open('rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
 status=dict(status='PARTIAL',input_samples=len(rows),input_barcodes=total,all_samples_LYPLA1_unique=True,all_sample_barcodes_unique=True,matrix_dimensions_and_entries_verified=True,
  source_sha256=h.hexdigest(),source_bytes=SOURCE.stat().st_size,malignancy_labels='NOT_YET_JOINED',differential_expression='NOT_RUN',note='GEO filtered CellRanger input; not the authors final139446-cell QC population. Mixed cell types; do not report as malignant LYPLA1 results.')
 (P/'preparation_validation.json').write_text(json.dumps(status,indent=2));print('PREPARATION_DONE',json.dumps(status),flush=True)
if __name__=='__main__':main()
