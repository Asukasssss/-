"""Read-only candidate suitability check, before spatial plotting or testing."""
import argparse,tarfile,io,h5py,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('archive');p.add_argument('--gene',default='LYPLA1');p.add_argument('--ensembl',default='ENSG00000120992');a=p.parse_args()
f=Path(a.archive)
with tarfile.open(f) as t:
 entries=[m for m in t.getmembers() if m.isfile() and Path(m.name).name=='filtered_feature_bc_matrix.h5'];assert len(entries)==1
 with h5py.File(io.BytesIO(t.extractfile(entries[0]).read()),'r') as h:
  names=[v.decode() for v in h['matrix/features/name'][:]];ids=[v.decode().split('.')[0] for v in h['matrix/features/id'][:]]
  out=dict(archive=f.name,gene=a.gene,ensembl=a.ensembl,symbol_matches=names.count(a.gene),ensembl_matches=ids.count(a.ensembl),features=len(names),filtered_barcodes=int(h['matrix/shape'][1]))
out['archive_sha256']=hashlib.sha256(f.read_bytes()).hexdigest();out['status']='PRESENT' if out['symbol_matches']==1 and out['ensembl_matches']==1 else 'NEEDS_REVIEW' if out['symbol_matches'] or out['ensembl_matches'] else 'NOT_EVALUABLE'
out['note']='Absent target is not zero expression; do not substitute LYPLA2 or LYPLAL1.'
print(json.dumps(out,indent=2))
