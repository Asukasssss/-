"""Independent aggregate checks for the actual SOP-adapted PDAC outputs."""
from pathlib import Path
import json,hashlib,subprocess
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
from prepare_sop_v3 import REPO,ROOT,RUNS,OLD_MET,OLD_PAT,OLD_MAP,TEMPLATES,dest,digest,write_json

def read(p):return pd.read_csv(p,sep='\t')
def main():
    checks={};allfamilies=[]
    oldmet=read(OLD_MET/'results.tsv')
    for name,mode in [('metabolite_paired.tsv','primary'),('metabolite_available_sensitivity.tsv','both_preimputation_available')]:
        d=read(dest('01_CAMP')/name).set_index('metabolite_key');o=oldmet[oldmet.analysis_type==mode].set_index('metabolite_key').loc[d.index]
        assert len(d)==307 and d.index.is_unique
        for c in ['effect','ci_lower','ci_upper','median_delta','rank_biserial']:
            np.testing.assert_allclose(d[c],o[c],equal_nan=True,rtol=0,atol=1e-12)
        np.testing.assert_allclose(d.p_original,o.p_value,equal_nan=True,rtol=0,atol=1e-15)
        ok=d.n>=8;np.testing.assert_allclose(d.loc[ok,'p_value'],o.loc[ok,'p_value'],rtol=0,atol=1e-15)
        assert d.loc[~ok,'p_value'].isna().all() and d.loc[~ok,'q_value'].isna().all()
        assert (d.n_up+d.n_down+d.n_equal==d.n).all()
        allfamilies.append(d)
    old=read(OLD_PAT/'associations.tsv')
    for name,mode,tier in [('tumor_association.tsv','primary','DIRECT'),('tumor_association_available.tsv','availability','DIRECT'),('conditional_tumor_association.tsv','primary','CONDITIONAL'),('conditional_tumor_association_available.tsv','availability','CONDITIONAL')]:
        d=read(dest('03_PATIENT')/name).set_index('relation_id');o=old[(old.analysis_type==mode)&(old.mapping_tier==tier)].set_index('relation_id').loc[d.index]
        assert d.index.is_unique
        for c in ['effect','p_value','n','seed','bootstrap_valid']:np.testing.assert_allclose(d[c],o[c],equal_nan=True,rtol=0,atol=1e-12)
        valid=d.bootstrap_valid>=3600
        for c in ['ci_lower','ci_upper']:np.testing.assert_allclose(d.loc[valid,c],o.loc[valid,c],equal_nan=True,rtol=0,atol=1e-12)
        allfamilies.append(d)
    for d in allfamilies:
        valid=d.p_value.notna();np.testing.assert_allclose(d.loc[valid,'q_value'],false_discovery_control(d.loc[valid,'p_value'].to_numpy()),rtol=0,atol=1e-12)
        assert d.loc[~valid,'q_value'].isna().all();assert (d.family_n_evaluable==valid.sum()).all()
    direct=read(dest('02_MAPPING')/'direct_relations.tsv');cond=read(dest('02_MAPPING')/'conditional_relations.tsv');pool=read(dest('02_MAPPING')/'gene_pool_current.tsv');history=read(dest('02_MAPPING')/'gene_pool_history_union.tsv')
    assert set(direct.stable_gene_id)==set(pool.stable_gene_id);assert len(pool)==250 and len(history)==687
    assert history.stable_gene_id.is_unique and history.gene.is_unique and set(pool.stable_gene_id)<=set(history.stable_gene_id)
    assert len(direct)==357 and len(cond)==358 and not pd.concat([direct,cond]).relation_id.duplicated().any()
    integrated=read(dest('07_INTEGRATION')/'candidate_relations_integrated.tsv');genes=read(dest('07_INTEGRATION')/'candidate_genes_integrated.tsv')
    assert len(integrated)==715 and integrated.relation_id.is_unique and len(genes)==687
    assert integrated[['RNA_p','RNA_q','RNA_effect']].isna().all().all();assert integrated.RNA_background.str.startswith('ACCESS_BLOCKED').all()
    assert len(read(dest('06_EXTERNAL')/'sc_gene_coverage.tsv'))==687*3
    templatechecks=0
    for stage in RUNS:
        for file in dest(stage).glob('*.tsv'):
            t=TEMPLATES/file.name
            if t.exists():
                prefix=list(read(t).columns);assert list(read(file).columns)[:len(prefix)]==prefix;templatechecks+=1
    # Verify all historic versioned files remain byte-identical in Git's normalized form.
    changes=subprocess.check_output(['git','diff','--name-only','109410a8f80537c428746e4551bdd760230a6e7f','--','results/PDAC'],cwd=REPO,text=True).splitlines()
    assert all(any(run in c for run in RUNS.values()) for c in changes),changes
    checks={'status':'PASS','six_BH_families_independently_checked_with_SciPy':True,'all_reused_effect_P_CI_checked':True,'metabolite_min8_gate_checked':True,'current_direct_gene_set_equality':250,'history_gene_union_unique':687,'relation_join_rows':715,'template_prefix_tables_checked':templatechecks,'server_numerics_executed':False,'matrix_input_current_hashes_reverified':False,'paired_t_and_SC_placeholders_not_misreported_as_results':True}
    for stage in RUNS:write_json(dest(stage)/'independent_validation.json',checks)
    print(json.dumps(checks))

if __name__=='__main__':main()
