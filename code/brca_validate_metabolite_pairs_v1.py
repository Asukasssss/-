"""Server validation of independent R signed-rank, counts and BH."""
from pathlib import Path
import argparse,json,hashlib
import pandas as pd,numpy as np
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);out=p.parse_args().out
d=pd.read_csv(out/'public/paired_metabolite318.tsv',sep='\t');r=pd.read_csv(out/'private/independent_R_checks.tsv',sep='\t');z=d.merge(r,on=['metabolite_name','test_family'],validate='one_to_one',suffixes=('','_R'))
assert len(z)==636 and z.p_value.notna().equals(z.p.notna())
err=dict(p=float(np.nanmax(abs(z.p_value-z.p))),q=float(np.nanmax(abs(z.q_value-z.q))))
assert max(err.values())<1e-10
assert z.n.eq(z.n_R).all() and z.pairs_higher.eq(z.higher).all() and z.pairs_lower.eq(z.lower).all() and z.pairs_equal.eq(z.equal).all()
for _,v in pd.read_csv(out/'public/source_manifest.tsv',sep='\t').iterrows():assert hashlib.sha256(Path(v.path).read_bytes()).hexdigest()==v.sha256
val=json.loads((out/'public/validation.json').read_text());val['independent_R_max_errors']=err;val['independent_R_direction_counts_match']=True
(out/'public/validation.json').write_text(json.dumps(val,indent=2))
pd.DataFrame([dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in [Path(__file__),out/'brca_check_metabolite_pairs_v1.R']]).to_csv(out/'public/verification_manifest.tsv',sep='\t',index=False)
print(json.dumps(err))
