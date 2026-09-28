"""Cell-weighted descriptive view, retaining the preceding CNA admission rules."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
PRIOR=ROOT/'results/collaborative/COAD/B/20260925T140257Z_uckl1_cna_v1'
sys.path.insert(0,str(PRIOR))
from run_uckl1_cna_v1 import build_labels,sha,PREFIX,GROUPS

ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
assert a.out.parent.resolve()==ROOT/'results/collaborative/COAD/B' and (a.out/'.running').is_dir()
pub=a.out/'public';pub.mkdir()
j,paths,coverage,audit=build_labels()
assert {k:sha(p) for k,p in paths.items()}==json.loads((PRIOR/'admission_hashes.json').read_text())
meta=pd.read_csv(paths['cell_metadata'],sep='\t',dtype=str,keep_default_na=False)
ix=pd.Index(meta.cell_id).get_indexer(j.cell_id)
assert (ix>=0).all() and len(set(ix))==len(meta)==len(j)
with np.load(paths['counts']) as z:counts=z['UCKL1'][ix];total=z['total'][ix]
cp10k=10000*counts/total;log=np.log1p(cp10k)
prior=pd.read_csv(PRIOR/'private_donor_values.tsv',sep='\t')
rows=[];hist=[]
for group in GROUPS:
    mask=(j.group==group).to_numpy();v=log[mask];raw=cp10k[mask];n=int(mask.sum())
    f=prior[prior.group==group]
    assert n==int(f.n_cells.sum())
    assert np.isclose(v.mean(),np.average(f.mean_log1p_CP10K,weights=f.n_cells),atol=1e-12)
    assert np.isclose((counts[mask]>0).mean(),np.average(f.detection_fraction,weights=f.n_cells),atol=1e-12)
    row=dict.fromkeys(PREFIX,'NA');row.update(cancer='COAD',cohort='Uhlitz_GSE166555',stage_id='06_EXTERNAL',run_id=a.out.name,
        analysis_version='COAD_UCKL1_cell_pooled_v1',analysis_type='cell_weighted_expression_description',gene='UCKL1',unit='cell',n=n,
        effect_type='cell_equal_mean_log1p_CP10K',effect=float(v.mean()),status='DONE',reason='DESCRIPTIVE_ONLY;CNN_NOT_PROVEN_NONMALIGNANT',
        source_id='Uhlitz author CNA calls;prior admission unchanged',group=group,n_donors=int(f.patient.nunique()),
        detected_cells=int((counts[mask]>0).sum()),detection_fraction=float((counts[mask]>0).mean()),
        mean_CP10K=float(raw.mean()),median_log1p_CP10K=float(np.median(v)),q75_log1p_CP10K=float(np.quantile(v,.75)),
        q90_log1p_CP10K=float(np.quantile(v,.90)),positive_cell_mean_log1p_CP10K=float(v[v>0].mean()))
    rows.append(row)
    edges=np.r_[0,np.arange(.1,10.1,.1),np.inf]
    h,_=np.histogram(v[v>0],bins=edges)
    hist.append(dict(group=group,bin='zero',lower=0,upper=0,n_cells=int((v==0).sum())))
    hist.extend(dict(group=group,bin=str(k),lower=float(edges[k]),upper=float(edges[k+1]),n_cells=int(c)) for k,c in enumerate(h))
    assert sum(x['n_cells'] for x in hist if x['group']==group)==n
pd.DataFrame(rows).to_csv(pub/'results.tsv',sep='\t',index=False)
pd.DataFrame(hist).to_csv(pub/'expression_histogram.tsv',sep='\t',index=False)
spec=dict(analysis_version='COAD_UCKL1_cell_pooled_v1',unit='cell',weight='each cell equal',gene='UCKL1',
    prior_run=PRIOR.name,prior_code_commit='0720928cafcbdcd864d606eecf873923e25424f9',
    admission='Exact prior groups and donor conflict exclusion;no new selection by expression',
    metrics=['mean log1p(CP10K)','mean CP10K','detection fraction','distribution quantiles'],
    tests='none;descriptive overall view requested by user',zero='retain measured zeros',
    interpretation='Cell-rich samples contribute more;not patient-equal or independent-cell inference')
(pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2)+'\n')
inputs=dict(paths,prior_donor_values=PRIOR/'private_donor_values.tsv',prior_script=PRIOR/'run_uckl1_cna_v1.py')
pd.DataFrame([dict(source_id=k,source_path=str(p),sha256=sha(p)) for k,p in inputs.items()]).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
(pub/'validation.json').write_text(json.dumps(dict(status='PASS',prior_input_hashes_identical=True,
    direct_cell_means_match_cell_count_weighted_donor_means=True,detection_counts_crosschecked=True,
    histogram_conserves_cells=True,no_new_p_or_q=True,all_groups_cell_median_zero=all(r['median_log1p_CP10K']==0 for r in rows)),indent=2)+'\n')
print(json.dumps(rows));(a.out/'DONE').write_text('DONE\n');(a.out/'.running').rmdir()
