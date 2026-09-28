"""Coverage-only gate: absent target never becomes a zero-expression observation."""
import argparse,hashlib,json
from pathlib import Path
import h5py,numpy as np
p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);a=p.parse_args();path=Path(a.input)
with h5py.File(path) as f:
 m=f['matrix'];names=m['features/name'][:].astype(str);ids=m['features/id'][:].astype(str)
 symbol=int(np.sum(names=='LYPLA1'));ensembl=int(np.sum(ids=='ENSG00000120992'));shape=m['shape'][:].tolist()
h=hashlib.sha256()
with path.open('rb') as f:
 for b in iter(lambda:f.read(1048576),b''):h.update(b)
result={'file':path.name,'sha256':h.hexdigest(),'bytes':path.stat().st_size,'gene':'LYPLA1','ensembl':'ENSG00000120992','symbol_matches':symbol,'ensembl_matches':ensembl,'matrix_shape':shape,'status':'NOT_EVALUABLE' if symbol==0 and ensembl==0 else 'NEEDS_REVIEW','scope':'Only this downloaded sample; no inference about other samples','reason':'Missing feature is not a measured zero; no substitution with paralogs; no expression test run'}
Path(a.output).write_text(json.dumps(result,indent=2));print(json.dumps(result))
