"""Server-only, locked 159-feature paired metabolite contrasts; public aggregates only."""
import argparse, csv, hashlib, json, platform
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from statsmodels.stats.multitest import multipletests
from patient_statistics import paired_rna, bh, stable_seed


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def save(p, x):
    with p.open('x', encoding='utf-8') as f:
        json.dump(x, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', type=Path, required=True)
    ap.add_argument('--run', type=Path, required=True)
    ap.add_argument('--commit', required=True)
    args = ap.parse_args()
    root, out = args.root.resolve(), args.run.resolve()
    assert out.parent == root / 'results/collaborative/COAD/B'
    assert (out / '.running').read_text().strip() == out.name
    spec = json.loads((out / 'locked_spec.json').read_text())
    assert spec['planned_features'] == 159 and spec['patients'] == 33
    res = root / 'results/collaborative/COAD/B/20260919T115842Z_identity_units_v1'
    clinpath = res / 'private_clinical_join.tsv'
    xp = root / 'data/candidates/camp_primary_tissue_multicancer/processed_metabolomics/PreprocessedData_COAD.xlsx'
    frozen = out / 'effects.tsv'
    assert sha(frozen) == '1bb1d57480a0fbc2185d11f7598e67e7443aef9e40ea8a036c6da4e5503a4597'
    inputs = [xp, clinpath, frozen, out/'locked_spec.json', out/'statistical_result_template.tsv', Path(__file__).resolve(), out/'patient_statistics.py']
    inputs += sorted(out.glob('*Metabo.R'))
    hashes = {str(p.relative_to(root)): sha(p) for p in inputs}
    save(out/'analysis_spec.json', dict(spec, run_id=out.name, code_commit=args.commit,
        written_before_contrasts_utc=datetime.now(timezone.utc).isoformat(), input_sha256=hashes,
        software=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__)))
    clin = pd.read_csv(clinpath, sep='\t', dtype=str, keep_default_na=False)
    t = clin[clin.TN.eq('Tumor') & clin.stage.isin(['stage I','stage II','stage III','stage IV'])]
    normals = clin[clin.TN.eq('Normal')]
    assert len(t)==33 and t.individual.is_unique and normals.individual.is_unique
    n = normals.set_index('individual').loc[t.individual]
    assert list(n.index)==list(t.individual)
    assert len(set(t.MetabID) | set(n.MetabID))==66
    sheets = pd.read_excel(xp, sheet_name=['data','data_imputed','metabo_imputed_filtered_Tumor','metabo_imputed_filtered_Normal'], index_col=0)
    for x in sheets.values():
        x.index=x.index.astype(str); x.columns=x.columns.astype(str)
        assert x.index.is_unique and x.columns.is_unique
    raw, full, mt, mn = [sheets[k] for k in ['data','data_imputed','metabo_imputed_filtered_Tumor','metabo_imputed_filtered_Normal']]
    assert full.loc[mt.index, mt.columns].equals(mt)
    assert full.loc[mn.index, mn.columns].equals(mn)
    observed = raw.stack()
    assert np.allclose(observed.values, full.stack().reindex(observed.index).values, rtol=0, atol=1e-10)
    with frozen.open(encoding='utf-8-sig', newline='') as f:
        effects=[r for r in csv.DictReader(f,delimiter='\t') if r['cancer']=='COAD']
    names=[r['feature_name'] for r in effects]
    assert len(names)==len(set(names))==159
    assert set(names)<=set(mt.index) and set(names)<=set(mn.index)
    prefix=(out/'statistical_result_template.tsv').read_text().strip().split('\t')
    extra=['original_effect','original_q','family_n_planned','up_pairs','down_pairs','equal_pairs','up_percent','down_percent',
           'both_observed_pairs','tumor_only_missing_pairs','normal_only_missing_pairs','both_missing_pairs',
           'n_nonzero_differences','bootstrap_valid','random_seed']
    records=[]; manual_p=[]; audited=0
    for original in effects:
        name=original['feature_name']
        tv=mt.loc[name,t.MetabID].to_numpy(float); nv=mn.loc[name,n.MetabID].to_numpy(float)
        rt=raw.loc[name,t.MetabID].to_numpy(float); rn=raw.loc[name,n.MetabID].to_numpy(float)
        assert np.isfinite(tv).all() and np.isfinite(nv).all()
        mask=np.isfinite(rt)&np.isfinite(rn)
        counts=dict(both_observed_pairs=int(mask.sum()),tumor_only_missing_pairs=int((~np.isfinite(rt)&np.isfinite(rn)).sum()),
            normal_only_missing_pairs=int((np.isfinite(rt)&~np.isfinite(rn)).sum()),both_missing_pairs=int((~np.isfinite(rt)&~np.isfinite(rn)).sum()))
        assert sum(counts.values())==33
        for kind, keep in [('primary33',np.ones(33,dtype=bool)),('observed_pairs',mask)]:
            x,y=tv[keep],nv[keep]; d=x-y
            result=paired_rna(x,y,stable_seed('paired_metabolites33_v1|'+name+'|'+kind),bootstrap=4000,minimum_n=8)
            row=dict.fromkeys(prefix+extra, 'NA')
            row.update(cancer='COAD',cohort='CAMP_COAD',stage_id='03_PATIENT' if kind=='primary33' else '04_ROBUSTNESS',
                run_id=out.name,analysis_version='paired_metabolites33_v1',analysis_type=kind,
                metabolite_key=original['metabolite_key'],metabolite_name=name,unit='GEO_annotated_patient_pair',
                effect_type='median_paired_T_minus_N_author_normalized_log_scale',test_family='COAD_paired_metabolites33_v1_'+kind,
                source_id=xp.name,original_effect=original['hedges_g'],original_q=original['effect_fdr'],family_n_planned=159,
                up_pairs=int((d>0).sum()),down_pairs=int((d<0).sum()),equal_pairs=int((d==0).sum()),
                up_percent=100*float((d>0).mean()) if len(d) else 'NA',down_percent=100*float((d<0).mean()) if len(d) else 'NA',**counts)
            row.update(result)
            assert row['up_pairs']+row['down_pairs']+row['equal_pairs']==row['n']
            # Independent direct comparisons and scalar median, without the subtraction masks.
            assert row['up_pairs']==sum(float(a)>float(b) for a,b in zip(x,y))
            assert row['down_pairs']==sum(float(a)<float(b) for a,b in zip(x,y))
            if result['status']=='DONE':
                import statistics
                assert result['effect']==statistics.median([float(a)-float(b) for a,b in zip(x,y)])
                # Independently reconstruct tie-corrected signed-rank normal approximation.
                z=d[d!=0]; ranks=stats.rankdata(abs(z)); m=len(z)
                if not m: p=1.0
                else:
                    _, ties=np.unique(abs(z),return_counts=True)
                    var=m*(m+1)*(2*m+1)/24 - sum(ties**3-ties)/48
                    p=float(2*stats.norm.sf(abs((ranks[z>0].sum()-m*(m+1)/4)/np.sqrt(var))))
                assert np.isclose(p,result['p_value'],rtol=1e-10,atol=1e-14)
                audited+=1
            records.append(row)
    for kind in ['primary33','observed_pairs']:
        family=[r for r in records if r['analysis_type']==kind]
        vals=[r['p_value'] if r['status']=='DONE' else 1.0 for r in family]
        qs=bh(vals)
        assert np.allclose(qs,multipletests(vals,method='fdr_bh')[1],rtol=1e-12,atol=1e-14)
        for row,q in zip(family,qs):
            row['family_n_evaluable']=sum(r['status']=='DONE' for r in family)
            if row['status']=='DONE': row['q_value']=float(q)
    for kind in ['primary33','observed_pairs']:
        with (out/(kind+'.tsv')).open('x',encoding='utf-8',newline='') as f:
            w=csv.DictWriter(f,prefix+extra,delimiter='\t',lineterminator='\n'); w.writeheader()
            w.writerows(r for r in records if r['analysis_type']==kind)
    summary={}
    for kind in ['primary33','observed_pairs']:
        rows=[r for r in records if r['analysis_type']==kind]
        sig=[r for r in rows if r['status']=='DONE' and r['q_value']<.05]
        summary[kind]=dict(planned=159,computable=sum(r['status']=='DONE' for r in rows),q_significant=len(sig),
            majority_up=sum(r['up_pairs']>r['n']/2 for r in sig),majority_down=sum(r['down_pairs']>r['n']/2 for r in sig),
            other_significant=sum(max(r['up_pairs'],r['down_pairs'])<=r['n']/2 for r in sig))
    summary['primary_pairs']=33
    summary['outside_frozen_family_T']=sorted(set(mt.index)-set(names))
    summary['outside_frozen_family_N']=sorted(set(mn.index)-set(names))
    summary['fully_observed_features']=sum(r['analysis_type']=='primary33' and r['both_observed_pairs']==33 for r in records)
    save(out/'summary.json',summary)
    with (out/'source_manifest.tsv').open('x',encoding='utf-8',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['source_path','version','sha256'])
        for p,h in hashes.items():w.writerow([p,'input_at_run',h])
    save(out/'validation.json',dict(status='PASS',frozen159_hash=True,pairs=33,unique_specimens=66,
        explicit_individual_join=True,T_N_equal_common_scale_subsets=True,observed_values_preserved=True,
        every_count_and_effect_independently_checked=True,wilcoxon_reconstructed_tests=audited,
        BH_statsmodels_check=True,output_sha256={p.name:sha(p) for p in [out/'primary33.tsv',out/'observed_pairs.tsv',out/'summary.json']},
        not_verified=['genotype-level patient identity','raw mass-spectrometry peak identity','logarithm base; no fold change reported']))
    save(out/'DONE.json',dict(status='DONE',run_id=out.name))
    (out/'.running').unlink()
    print(json.dumps(summary))


if __name__=='__main__':
    main()
