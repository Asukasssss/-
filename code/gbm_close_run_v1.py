from pathlib import Path
import json,hashlib,pandas as pd,datetime
r=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/GBM/A/20260925T103000Z_discovery_v1')
p=r/'public/07_INTEGRATION';d=pd.read_csv(p/'checksums.tsv',sep='\t')
for _,z in d.iterrows():assert hashlib.sha256((p/str(z['file']).replace(chr(92),'/')).read_bytes()).hexdigest()==z.sha256
assert (r/'.running').read_text()=='gbm_discovery_v1'
summary=dict(execution_status='COMPLETED_CURRENT_SCOPED_RUN',scientific_status='PARTIAL_MAPPING457_PENDING',UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),server_integration_checksums=len(d),source_patient_matrices_retained_on_server=True)
(r/'execution_summary.json').write_text(json.dumps(summary,indent=2));(r/'.running').rename(r/'.completed_scoped');print(json.dumps(summary))
