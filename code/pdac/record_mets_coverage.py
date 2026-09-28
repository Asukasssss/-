"""Record the observed Spatial-assay gate and the later inspection failure separately."""
import argparse,json,hashlib,re,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--run',required=True);p.add_argument('--source',required=True);p.add_argument('--download-log',required=True);a=p.parse_args();run=Path(a.run);src=Path(a.source)
log=(run/'inspect.log').read_text();download=Path(a.download_log).read_text()
m=re.search(r'ASSAY Spatial Assay (\d+) (\d+)\s*\nTARGET[ \t]*\r?\n',log)
assert m,'No explicit empty Spatial target result'
assert 'Seurat' in log and 'Error' in log,'Expected subsequent dependency limitation must remain recorded'
assert '|OK' in download.replace(' ','') and not Path(str(src)+'.aria2').exists(),'Download not complete'
h=hashlib.sha256()
with src.open('rb') as f:
 for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
result={'status':'NOT_EVALUABLE','cohort':'GSE272362','gene':'LYPLA1','ensembl':'ENSG00000120992','assay':'Spatial','assay_features':int(m[1]),'assay_spots':int(m[2]),'symbol_or_ensembl_matches_in_assay_rownames':0,'scope':'Author Spatial assay rownames only; original count slot independently checking NOT_PERFORMED; other assay layers not inspected successfully','inspection_status':'PARTIAL: later SCT assay dispatch requires unavailable Seurat package','expression_tests':'NOT_RUN','reason':'Published Spatial assay rownames have no target match. This is not zero expression or biological absence; do not extrapolate to uninspected original sources.','source':{'url':'https://zenodo.org/records/22695225/files/PDAC_Updated_V2.rds','version':'22695225 v2','bytes':src.stat().st_size,'sha256':h.hexdigest(),'download_expected_md5':'024d9ce0dd48abf3bd425c7ebead03fc','aria2_checksum_verified_download_status':'OK'}}
pub=run/'public';(pub/'target_coverage.json').write_text(json.dumps(result,indent=2));shutil.copyfile(run/'analysis_spec.json',pub/'analysis_spec.json');shutil.copyfile(run/'inspect.log',pub/'inspection_log.txt')
(pub/'validation.json').write_text(json.dumps({'source_download_checksum':'PASS','Spatial_assay_target_presence':'ABSENT_IN_ASSAY_ROWNAMES','subsequent_object_inspection':'FAILED_MISSING_Seurat_PACKAGE','raw_counts_independent_gate':'NOT_PERFORMED','malignant_enrichment':'NOT_EVALUABLE','no_expression_recalculation':True},indent=2))
(run/'.running').rmdir();print(json.dumps(result))
