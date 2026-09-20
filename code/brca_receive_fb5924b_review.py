"""Receive aggregate review; compare excerpts to pinned Git objects, no patient tests."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import zipfile
import numpy as np
import pandas as pd

BASE = 'results/BRCA/04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2/'

def frozen(name):
    return pd.read_csv(io.BytesIO(subprocess.check_output(['git', 'show', 'fb5924b:'+BASE+name])), sep='\t')

def main(archive, out):
    out.mkdir(parents=True, exist_ok=False)
    z = zipfile.ZipFile(archive)
    prefix = z.namelist()[0].split('/')[0]+'/'
    manifest = pd.read_csv(io.BytesIO(z.read(prefix+'PACKAGE_SHA256.tsv')), sep='\t')
    for r in manifest.itertuples():
        data = z.read(prefix+r.path)
        assert len(data)==r.bytes and hashlib.sha256(data).hexdigest()==r.sha256, r.path
    # Only named aggregate evidence; never execute supplied code automatically.
    for name in ['00_README_CN.md','PACKAGE_SHA256.tsv','validation.json',
                 'code/verify_review.py','sources/source_registry.tsv',
                 'tables/association_excerpts_14.tsv','tables/source_stability_excerpts_17.tsv',
                 'tables/review_findings.tsv']:
        dest=out/'imported_review'/name
        dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_bytes(z.read(prefix+name))
    assoc=pd.read_csv(out/'imported_review/tables/association_excerpts_14.tsv',sep='\t')
    source=pd.read_csv(out/'imported_review/tables/source_stability_excerpts_17.tsv',sep='\t')
    truth=frozen('all174_composition_sensitivity.tsv').set_index(['relation_id','test_family'],verify_integrity=True)
    sc=frozen('all117_next_actions_CN.tsv').set_index('gene',verify_integrity=True)
    checks=[]
    def check(key,field,actual,expected):
        ok=bool(np.isclose(float(actual),float(expected),rtol=0,atol=1e-12)) if isinstance(expected,(int,float,np.number)) else str(actual)==str(expected)
        checks.append(dict(key=key,field=field,review_value=actual,pinned_value=expected,status='DONE' if ok else 'NEEDS_REVIEW'))
    for _,r in assoc.iterrows():
        for model in ['ER_PC12','ER_MonoEndoFib']:
            t=truth.loc[(r.relation_id,model)]
            for a,b in {'rho':'effect','ci_lower':'ci_lower','ci_upper':'ci_upper','p':'p_value','q':'q_value'}.items():
                check(r.relation_id+'|'+model,a,r[model+'_'+a],t[b])
            for a,b in {'n_specimens':'n','ER_same_subset_rho':'ER_v3_rho','ER_same_subset_q':'ER_v3_q','gene':'gene'}.items():
                check(r.relation_id+'|'+model,a,r[a],t[b])
    for _,r in source.iterrows():
        t=sc.loc[r.gene]
        for short,cohort in [('Wu','Wu2021'),('Pal','Pal2021_reprocessed')]:
            for a,b in {'top':'top_lineage','bootstrap_top_frequency':'bootstrap_top_frequency','runner':'runner_lineage','paired_n':'paired_source_labels','paired_positive_fraction':'paired_positive_fraction'}.items():
                check(r.gene,short+'_'+a,r[short+'_'+a],t['rob_'+cohort+'_'+b])
        check(r.gene,'two_source_stable_rule',r.two_source_stable_rule,t.rob_two_tumor_source_stable_descriptive)
    pd.DataFrame(checks).to_csv(out/'excerpt_checks.tsv',sep='\t',index=False,lineterminator='\n')
    receipt={'status':'DONE' if all(x['status']=='DONE' for x in checks) else 'NEEDS_REVIEW',
             'source_commit':subprocess.check_output(['git','rev-parse','fb5924b'],text=True).strip(),
             'archive_sha256':hashlib.sha256(archive.read_bytes()).hexdigest(),
             'package_hashes_verified':len(manifest),'association_relations':len(assoc),
             'source_genes':len(source),'scalar_checks':len(checks),
             'mismatches':sum(x['status']!='DONE' for x in checks),
             'patient_tests':0,'full_BH_rerun':False,'full_literature_reaudit':False,
             'source_citation_note':'Use all174_composition_sensitivity.tsv for relation-level comparisons; best-gene summaries may refer to different metabolites.'}
    (out/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--archive',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();main(a.archive,a.out)
