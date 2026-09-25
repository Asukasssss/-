"""Hash server-resident input files actually used, without exporting matrices."""
from pathlib import Path
import hashlib,json
import pandas as pd
R=Path(__file__).resolve().parent
D=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates')
idx=pd.read_csv(R/'private/units_index.tsv',sep='\t',dtype={'unit':str})
paths=[D/'gse263733_pdac_scRNA_v0.1/GSE263733_Raw_counts.txt.gz']
paths+=list((D/'gse278688_pdac_scRNA_v0.1').glob('GSE278688_sc_*.gz'))
for unit in idx.loc[idx.cohort=='GSE242230','unit']:
    a=list((D/'GSE242230/raw_matrices').glob('*_'+unit+'_matrix.mtx.gz'));assert len(a)==1
    prefix=str(a[0])[:-len('matrix.mtx.gz')]
    paths.extend([a[0],Path(prefix+'barcodes.tsv.gz'),Path(prefix+'features.tsv.gz')])
paths+=list(R.glob('myeloid_infercnv_*'))+[R/'public/analysis_spec.json',R/'private/gene_order_no_chr3.tsv']
rows=[]
for p in sorted(set(paths)):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    rows.append({'path':str(p),'sha256':h.hexdigest()})
pd.concat([pd.read_csv(R/'public/source_manifest.tsv',sep='\t'),pd.DataFrame(rows)]).drop_duplicates('path',keep='last').to_csv(R/'public/source_manifest.tsv',sep='\t',index=False)
(R/'INPUT_HASHES_DONE.json').write_text(json.dumps({'new_files_hashed':len(rows)}))
print('INPUT_HASHES_DONE',len(rows),flush=True)
