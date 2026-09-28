"""Audit exclusion of LYPLA1 in the older probe-based official HD showcase."""
import h5py,requests,json,hashlib
from brca_hd_range_v1 import RemoteH5,URL
from brca_gse210616_stream_to_server import remote,transfer
from brca_hd_fetch_v1 import D,OUT
r=RemoteH5(URL)
with h5py.File(r,'r') as f:
 g=f['features'];i=list(g['name'][:]).index(b'LYPLA1');sets={k:bool(i in g['target_sets'][k][:]) for k in g['target_sets']}
 result={'dataset':'10x Visium HD FF Human Breast Cancer (probe-based DCIS)','official_feature_slice':URL,'feature_index':i,'symbol':g['name'][i].decode(),'ensembl':g['id'][i].decode(),'included_in_target_sets':sets,'raw_feature_positions':int(f['feature_slices'][str(i)]['data'].shape[0]),'interpretation':'Raw counts exist but gene excluded from retained target set. Do not use raw counts for reliable LYPLA1 spatial expression.','http_range_bytes':r.transferred}
u='https://raw.githubusercontent.com/ejscience/2025_Wamaitha_Rhesus_OvarianReserve/main/Visium_SpaceRanger/RhesusMacaqueProbeLists/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv'
p=requests.get(u,timeout=30);p.raise_for_status();lines=[s for s in p.text.splitlines() if 'ENSG00000120992,' in s];assert len(lines)==1 and ',FALSE,' in lines[0]
result.update(probe_reference_mirror=u,probe_reference_sha256=hashlib.sha256(p.content).hexdigest(),LYPLA1_probe_included=False,official_exclusion_documentation='https://www.10xgenomics.com/support/spatial-gene-and-protein-expression/documentation/steps/probe-sets/visium-human-transcriptome-probe-set-v2-0')
transfer([p.content],D+'/probe_reference_v2.csv');transfer([json.dumps(result,indent=2).encode()],OUT+'/public/probe_coverage_audit.json');print(json.dumps(result,indent=2))
