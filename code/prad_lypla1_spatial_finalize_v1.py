"""Standard-format cohort summaries and completion checks; run after analysis.
Usage: python prad_lypla1_spatial_finalize_v1.py RUN_ROOT BASE_GIT_COMMIT
Does not rerun statistical tests or release individual section measurements.
"""
from pathlib import Path
import sys,json,hashlib
import pandas as pd
root=Path(sys.argv[1]);pub=root/'public'
c=pd.read_csv(root/'private/section_contrasts.tsv',sep='\t')
g=pd.read_csv(pub/'tumor_celltype_summary.tsv',sep='\t')
assert g.n_beads.sum()==83632 and len(c)==12
assert c.groupby('comparison').size().eq(4).all()
for grid in [4,8]:
    assert (abs(c[f'rate_ratio_grid{grid}']-c.pseudobulk_ratio)<1e-8).all()
records=[]
for comp,x in c.groupby('comparison'):
    records.append(dict(cancer='PRAD',cohort='Hirz2023_SlideSeqV2',stage_id='06_EXTERNAL',run_id=root.name,
        analysis_version='prad_lypla1_spatial_v1',analysis_type='spatial_author_label_enrichment',metabolite_key='NA',metabolite_name='NA',gene='LYPLA1',
        unit='section_summary_from_two_cancer_donors',n=4,n_reference=2,effect_type='section_equal_mean_depth_offset_rate_ratio',effect=x.pseudobulk_ratio.mean(),
        ci_lower='NA',ci_upper='NA',p_value='NA',q_value='NA',test_family='12 bead rank tests; 24 spatial sensitivity tests',family_n_evaluable=36,
        status='DONE',reason='Descriptive section-equal effect; no population P estimated with two donors; per-section nominal P summarized separately',
        source_id='doi:10.1038/s41467-023-36325-2',comparison=comp,
        min_section_effect=x.pseudobulk_ratio.min(),max_section_effect=x.pseudobulk_ratio.max()))
pd.DataFrame(records).to_csv(pub/'results.tsv',sep='\t',index=False)
v=json.loads((pub/'validation.json').read_text());v.update(status='DONE',evidence_level='EXPLORATORY',
    glm_rate_equals_pseudobulk_ratio=True,public_group_sum_verified=True,base_git_commit=sys.argv[2])
(pub/'validation.json').write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n')
s=json.loads((pub/'analysis_spec.json').read_text());s['base_git_commit']=sys.argv[2]
s['code_sha256']={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in root.glob('prad_lypla1_spatial*') if f.is_file()}
(pub/'analysis_spec.json').write_text(json.dumps(s,ensure_ascii=False,indent=2)+'\n')
print('Final checks passed; public output contains cohort summaries and figures only.')
