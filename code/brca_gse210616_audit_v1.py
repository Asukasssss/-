"""Audit completeness, explicit patient mapping, file integrity and figure coverage."""
from pathlib import Path
import json,hashlib,re,sys,gzip
import pandas as pd,numpy as np
from PIL import Image
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/BRCA_GSE210616')
out=Path(sys.argv[1]);v=json.loads((out/'validation.json').read_text());records=json.loads((ROOT/'manifest.json').read_text())
assert len(records)==43 and v['sections_analyzed']==43 and not v['failures']
sec=pd.read_csv(out/'section_results.tsv',sep='\t');pat=pd.read_csv(out/'patient_results.tsv',sep='\t')
assert len(sec)==43 and sec.gsm.is_unique and len(pat)==22
mapped=[];bounds=[]
for r in records:
    match=re.search(r'patient (\d+), section (\d+)',r['title']);patient,section=map(int,match.groups());mapped.append((patient,section))
    rr=sec.loc[sec.gsm==r['gsm']].iloc[0];assert int(rr.patient)==patient and int(rr.section)==section
    for f in r['files']:assert hashlib.sha256(Path(f['path']).read_bytes()).hexdigest()==f['sha256']
    d=ROOT/r['gsm'];scale=json.load(gzip.open(next(d.glob('*scalefactors*.gz')),'rt'))['tissue_hires_scalef']
    pos=pd.read_csv(next(d.glob('*positions*.gz')),header=None).set_index(0)
    spots=pd.read_csv(out/f"{r['gsm']}_spot_measurements.tsv.gz",sep='\t');p=pos.loc[spots.barcode]
    assert len(spots)==int(rr.n_qc) and spots.barcode.is_unique
    im=Image.open(gzip.open(next(d.glob('*image*.gz')),'rb'));w,h=im.size
    x=p[5].to_numpy()*scale;y=p[4].to_numpy()*scale
    assert ((x>=0)&(x<w)&(y>=0)&(y<h)).all()
    assert (out/'figures'/f"{r['gsm']}.png").exists()
assert len(set(mapped))==43 and len(set(x[0] for x in mapped))==22
check={'status':'PASS','sections':43,'patients':22,'input_files_hashed':sum(len(r['files']) for r in records),
       'explicit_GEO_patient_section_identity':'PASS','all_QC_spots_inside_image_bounds':'PASS','all_43_figures_present':'PASS',
       'matrix_data_local_storage':False,'pathology_labels':'unavailable; not inferred as malignant','treatment':'unavailable; not imputed',
       'analysis_code_commit':'e75a07a','analysis_script_sha256':hashlib.sha256((out/'analyze.py').read_bytes()).hexdigest(),
       'execution_environment':'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1'}
(out/'audit.json').write_text(json.dumps(check,indent=2));print(json.dumps(check))
