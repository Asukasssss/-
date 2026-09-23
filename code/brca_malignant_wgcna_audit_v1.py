"""Server-side structural audit; publishes counts/hashes, never donor/cell rows."""
from pathlib import Path
import sys,json,hashlib,platform
import pandas as pd
import numpy as np
import scipy,h5py,sklearn
R=Path(sys.argv[1]);P=R/'public';D=R/'private'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
checks=[];manifest=[]
for co in ['Wu2021','Pal2021_reprocessed']:
    x=pd.read_csv(D/(co+'_expr.tsv'),sep='\t',index_col=0)
    meta=pd.read_csv(D/(co+'_metacells.tsv'),sep='\t')
    cell=pd.read_csv(D/(co+'_cell_membership.tsv'),sep='\t')
    donor=pd.read_csv(D/(co+'_donor_expr.tsv'),sep='\t',index_col=0)
    assert list(x.index)==meta.metacell.tolist()
    assert not cell.barcode.duplicated().any()
    assert meta.groupby('donor').size().eq(8).all()
    assert cell.groupby('metacell').size().eq(15).all()
    assert cell.groupby('metacell').donor.nunique().eq(1).all()
    assert set(cell.metacell)==set(meta.metacell)
    assert np.isfinite(x.to_numpy()).all()
    assert x.columns.is_unique and donor.columns.is_unique
    assert set(donor.index)==set(meta.donor)
    checks.append(dict(cohort=co,donors=len(donor),metacells=len(meta),genes=x.shape[1],cells=len(cell),zero_overlap=True,exactly_15_cells_per_metacell=True,eight_metacells_per_donor=True,unique_donor_per_metacell=True,finite_expression=True))
for p in sorted(list(R.glob('*.py'))+list(R.glob('*.R'))+list(R.glob('*config.json'))):manifest.append(dict(kind='executed_code_or_config',file=p.name,sha256=sha(p)))
for p in sorted(D.glob('*.tsv')):manifest.append(dict(kind='private_input_or_intermediate_NOT_EXPORTED',file=p.name,sha256=sha(p)))
pd.DataFrame(manifest).to_csv(P/'execution_manifest.tsv',sep='\t',index=False)
(P/'server_audit.json').write_text(json.dumps(dict(status='PASS',checks=checks,source_identity='publisher labels;no new genotype audit',software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,sklearn=sklearn.__version__,h5py=h5py.__version__)),indent=2))
print(json.dumps(checks))
