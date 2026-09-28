"""Check one complete sample independently using scipy MatrixMarket reader."""
from pathlib import Path
import tarfile,gzip,json
import pandas as pd,numpy as np
from scipy.io import mmread
R=Path(__file__).resolve().parent
src='/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/GSE205013/GSE205013_RAW.tar'
with tarfile.open(src) as t:
 name=next(m.name for m in t.getmembers() if m.name.endswith('_P04_matrix.mtx.gz'));prefix=name[:-len('_matrix.mtx.gz')]
 x=mmread(gzip.GzipFile(fileobj=t.extractfile(name))).tocsr()
 genes=pd.read_csv(gzip.GzipFile(fileobj=t.extractfile(prefix+'_features.tsv.gz')),sep='\t',header=None)
 bc=pd.read_csv(gzip.GzipFile(fileobj=t.extractfile(prefix+'_barcodes.tsv.gz')),sep='\t',header=None)[0]
d=pd.read_csv(R/'private/P04_LYPLA1.tsv.gz',sep='\t');ix=np.flatnonzero(genes[1].eq('LYPLA1'))
assert len(ix)==1 and list(d.barcode)==list(bc)
assert np.array_equal(x[ix[0]].toarray().ravel(),d.LYPLA1_count)
assert np.array_equal(np.asarray(x.sum(axis=0)).ravel(),d.total_UMI)
assert np.allclose(np.log1p(10000*d.LYPLA1_count/d.total_UMI),d.LYPLA1_log1p)
v=dict(status='PASS',scope='One complete sample independently read via scipy.io.mmread',n_cells=len(d),barcode_order_match=True,all_LYPLA1_counts_match=True,all_library_sizes_match=True,all_log1p_values_match=True,not_validated='Malignancy labels, final author QC, and clinical contrasts')
(R/'public/preparation_independent_check.json').write_text(json.dumps(v,indent=2));print(json.dumps(v))
