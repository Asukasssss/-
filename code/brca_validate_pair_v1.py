"""Validate completed server sensitivity outputs against independent R and provenance."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import pandas as pd

p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();out=a.out
assert (out/'NUMERICAL_DONE').exists()
e=pd.read_csv(out/'public/paired_RNA117.tsv',sep='\t').set_index('gene')
r=pd.read_csv(out/'private/R_paired_check.tsv',sep='\t').set_index('gene').loc[e.index]
assert e.p_value.notna().equals(r.p.notna())
errors={}
for x,y in [('effect','effect'),('p_value','p'),('q_value','q'),('ci_lower','lower'),('ci_upper','upper')]:
    errors[x]=float(np.nanmax(abs(e[x]-r[y])))
    # R qt and this SciPy t.ppf differ by 4.53e-9 at df44, yielding <=1.93e-9 CI differences.
    # Effects and P/q retain the original stricter tolerance; statistical results are unchanged.
    assert errors[x] < (1e-8 if x in ['ci_lower','ci_upper'] else 1e-10)
manifest=pd.read_csv(out/'public/source_manifest.tsv',sep='\t')
for _,v in manifest.iterrows():assert hashlib.sha256(Path(v.path).read_bytes()).hexdigest()==v.sha256
identity=out.parent/'20260921T102429Z_camp_sample_identity_v1/public/source_manifest.tsv'
old=pd.read_csv(identity,sep='\t');old=old[old.path.str.contains('transcriptomics_processed|processed_metabolomics|MasterMapping')]
assert len(old)==3
for _,v in old.iterrows():assert hashlib.sha256(Path(v.path).read_bytes()).hexdigest()==v.sha256
val=json.loads((out/'public/validation.json').read_text());val['independent_R_paired_max_errors']=errors;val['identity_audit_source_hashes_match']=3
val['R_check_interval_tolerance_note']='R/SciPy inverse t implementations differ at df44; CI absolute tolerance1e-8, effect/P/q1e-10; no estimates changed'
(out/'public/validation.json').write_text(json.dumps(val,indent=2))
files=[Path(__file__),out/'brca_check_paired_RNA_v1.R']
pd.DataFrame([dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files]).to_csv(out/'public/verification_manifest.tsv',sep='\t',index=False)
print(json.dumps(errors))
