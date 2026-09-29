"""Public provenance only; source data remain on server165."""
from pathlib import Path
import sys,hashlib,json,platform
import pandas as pd,numpy,h5py,scipy,matplotlib
R=Path(sys.argv[1]);B=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716');C=B/'results/collaborative'
files=[C/'BRCA/A/20260919T140105Z_scRNA117_v1/source/wu_curated.h5ad',B/'data/candidates/SCP1089/GSE202051_totaldata-final-toshare.h5ad',C/'PRAD/B/20260922T142000Z_discovery_v1/source/PRAD24_cellxgene.h5ad']
files += [C/'COAD/B/20260922T141755Z_source_contract_sc_v2'/n for n in ['private_cell_counts.npz','private_cell_metadata.tsv']]
files += [C/'COAD/B/20260925T140257Z_uckl1_cna_v1/public/results.tsv',C/'PRAD/B/20260925T133653Z_slc6a6_sc_paired_v1/public/06_EXTERNAL/results.tsv',C/'PRAD/B/20260925T133653Z_slc6a6_sc_paired_v1/private/paired_donor_values_private.tsv']
files += [C/'BRCA/A/20260929T071500Z_four_sc_v3'/(c+s) for c in ['BRCA','COAD','PDAC','PRAD'] for s in ['_private.npz','_meta.json']]
files += [B/'data/candidates/PDAC_GeoMx_Bell2025'/n for n in ['metadata_with_VI_subtypes_and_NGS.csv','ProbeQC_merged_batches.csv','fully_batch_corrected_vsd.csv']]
rows=[]
for p in files:
 assert p.exists(),p
 sha='NOT_REHASHED_LARGE_SOURCE'
 if p.stat().st_size<100000000:
  sha=hashlib.sha256(p.read_bytes()).hexdigest()
 rows.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha,version='existing source; read only'))
pd.DataFrame(rows).to_csv(R/'public/source_manifest.tsv',sep='\t',index=False)
versions=dict(python=platform.python_version(),pandas=pd.__version__,numpy=numpy.__version__,h5py=h5py.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__)
(R/'public/server_versions.json').write_text(json.dumps(versions,indent=2))
print('PROVENANCE_READY',len(rows))
