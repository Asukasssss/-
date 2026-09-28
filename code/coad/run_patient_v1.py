"""Run locked COAD patient families on server165; export only aggregate statistics."""
import argparse,csv,hashlib,json,platform,sys,traceback
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from patient_statistics import association,paired_rna,bh,stable_seed

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,obj):
    with p.open('x',encoding='utf-8') as f:json.dump(obj,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def table(p,fields,rows):
    with p.open('x',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for row in rows:w.writerow({k:'NA' if row.get(k) is None else row.get(k,'NA') for k in fields})

def main():
    p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,required=True);p.add_argument('--run-dir',type=Path,required=True);p.add_argument('--resolution-dir',type=Path,required=True);p.add_argument('--code-commit',required=True);a=p.parse_args()
    root,out,res=a.project_root.resolve(),a.run_dir.resolve(),a.resolution_dir.resolve()
    assert str(root).startswith('/public3/') and out.parent==res.parent==root/'results/collaborative/COAD/B'
    assert (out/'.running').read_text().strip()==out.name and not (out/'summary.json').exists()
    spec=json.loads((out/'locked_spec.json').read_text())
    cohort=json.loads((res/'cohort_summary.json').read_text());assert cohort['unit_check']=='DONE'
    clin=pd.read_csv(res/'private_clinical_join.tsv',sep='\t',dtype=str,keep_default_na=False)
    tm=clin[clin.TN.eq('Tumor')].copy();nm=clin[clin.TN.eq('Normal')].copy()
    assert len(tm)==37 and tm.individual.is_unique and nm.individual.is_unique
    primary=tm[tm.stage.isin(['stage I','stage II','stage III','stage IV'])].copy();assert len(primary)==33
    paired=nm.set_index('individual').loc[primary.individual]
    assert list(paired.index)==list(primary.individual)
    base=root/'data/candidates/camp_primary_tissue_multicancer'
    mp=base/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv'
    mapping=pd.read_csv(mp,dtype=str);mapping=mapping[mapping.Dataset.eq('COAD')]
    assert set(tm.RNAID)==set(mapping[mapping.TN.eq('Tumor')].RNAID)
    rp=base/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/mapping.RNAFile.iloc[0]
    xp=base/'processed_metabolomics'/mapping.MetabFile.iloc[0]
    inputs=[mp,rp,xp,res/'private_clinical_join.tsv',res/'candidate_eligibility.tsv',res/'cohort_summary.json',res/'geo_samples.soft',out/'effects.tsv',out/'locked_spec.json',out/'statistical_result_template.tsv',Path(__file__),out/'patient_statistics.py']
    hashes={p:sha(p) for p in inputs}
    # Check the public family and clinical join against the earlier unit-resolution run.
    prior=json.loads((res/'validation.json').read_text())
    assert sha(res/'candidate_eligibility.tsv')==prior['public_output_sha256']['candidate_eligibility.tsv']
    assert sha(res/'cohort_summary.json')==prior['public_output_sha256']['cohort_summary.json']
    rna=pd.read_csv(rp,index_col=0)
    sheet=mapping[mapping.TN.eq('Tumor')].MetabFile_sheet.unique();assert len(sheet)==1
    matrices=pd.read_excel(xp,sheet_name=[sheet[0],'data'],index_col=0)
    met,observed=matrices[sheet[0]],matrices['data']
    for matrix in [rna,met,observed]:
        matrix.columns=matrix.columns.astype(str);matrix.index=matrix.index.astype(str)
        assert matrix.index.is_unique and matrix.columns.is_unique
    assert tm.MetabID.isin(met.columns).all() and tm.RNAID.isin(rna.columns).all()
    assert primary.MetabID.isin(observed.columns).all() and paired.RNAID.isin(rna.columns).all()
    effects={(r['feature_name'],r['metabolite_key']):r for r in read(out/'effects.tsv') if r['cancer']=='COAD'}
    candidates=read(res/'candidate_eligibility.tsv')
    selected=[r for r in candidates if r['in_prespecified_family']=='True']
    assert len(candidates)==spec['all_candidate_pairs']==974 and len(selected)==spec['association_family_planned']==674
    genes={r['human_gene_id']:r for r in selected};assert len(genes)==spec['rna_family_planned']==458
    prefix=(out/'statistical_result_template.tsv').read_text().strip().split('\t')
    extra=['human_gene_id','rna_symbol','original_effect','original_q','family_n_planned','biochemical_disposition',
           'n_unique_metabolite','n_unique_rna','bootstrap_valid','permutations','random_seed','loo_valid','loo_min','loo_max','loo_max_abs_delta','loo_any_sign_change','n_nonzero_differences']
    fields=prefix+extra
    def base_row(r,kind):
        effect=effects.get((r['metabolite_name'],r['metabolite_key']))
        row=dict.fromkeys(fields,None)
        row.update(cancer='COAD',cohort='COAD',stage_id='04_ROBUSTNESS' if kind in ['all37','availability33'] else '03_PATIENT',
          run_id=out.name,analysis_version=spec['version'],analysis_type=kind,metabolite_key=r['metabolite_key'],metabolite_name=r['metabolite_name'],
          gene=r['gene'],human_gene_id=r['human_gene_id'],rna_symbol=r['rna_symbol'],unit='GEO_annotated_independent_individual',
          effect_type='Spearman_rho',test_family='COAD_patient_v1_'+kind,
          original_effect=effect['hedges_g'] if effect else None,original_q=effect['effect_fdr'] if effect else None,
          source_id=rp.name+';'+xp.name,biochemical_disposition=r.get('biochemical_disposition'),family_n_planned=674)
        return row
    records={k:[] for k in ['primary33','all37','availability33','paired_RNA33']}
    spec_out=dict(spec,created_utc=datetime.now(timezone.utc).isoformat(),code_commit=a.code_commit,
       run_id=out.name,input_hashes={str(p.relative_to(root)):h for p,h in hashes.items()},
       software={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__})
    dump(out/'analysis_spec.json',spec_out)  # Written before any association or contrast.
    for i,r in enumerate(candidates):
        for kind,people in [('primary33',primary),('all37',tm),('availability33',primary)]:
            row=base_row(r,kind)
            if r['in_prespecified_family']!='True':
                row.update(status='NEEDS_REVIEW',reason='OUTSIDE_LOCKED_FAMILY_'+r['biochemical_disposition'],test_family='NOT_IN_CURRENT_FAMILY',family_n_planned=None)
            elif r['rna_symbol']=='NA':row.update(status='NOT_EVALUABLE',reason='RNA_SYMBOL_UNAVAILABLE')
            else:
                x=pd.to_numeric(met.loc[r['metabolite_name'],people.MetabID],errors='coerce').to_numpy(dtype=float)
                y=pd.to_numeric(rna.loc[r['rna_symbol'],people.RNAID],errors='coerce').to_numpy(dtype=float)
                if kind=='availability33':
                    mask=np.isfinite(pd.to_numeric(observed.loc[r['metabolite_name'],people.MetabID],errors='coerce').to_numpy(dtype=float))
                    x=np.where(mask,x,np.nan)
                seed=stable_seed('|'.join([spec['version'],kind,r['metabolite_key'],r['human_gene_id']]))
                row.update(association(x,y,seed,permutations=spec['permutations'],bootstrap=spec['bootstrap_replicates'],minimum_n=spec['minimum_n'],loo=kind=='primary33'))
            records[kind].append(row)
        if (i+1)%50==0:print('RELATIONS',i+1,'/',len(candidates),flush=True)
    for ident,r in sorted(genes.items()):
        dummy=dict(r,metabolite_name='NA',metabolite_key='NA')
        row=base_row(dummy,'paired_RNA33');row.update(effect_type='median_paired_difference_author_RNA_scale',family_n_planned=458,unit='explicit_GEO_individual_matched_pair')
        if r['rna_symbol']=='NA':row.update(status='NOT_EVALUABLE',reason='RNA_SYMBOL_UNAVAILABLE')
        else:
            t=pd.to_numeric(rna.loc[r['rna_symbol'],primary.RNAID],errors='coerce').to_numpy(dtype=float)
            n=pd.to_numeric(rna.loc[r['rna_symbol'],paired.RNAID],errors='coerce').to_numpy(dtype=float)
            row.update(paired_rna(t,n,stable_seed('|'.join([spec['version'],'paired_RNA33','NA',ident])),bootstrap=spec['bootstrap_replicates'],minimum_n=spec['minimum_n']))
        records['paired_RNA33'].append(row)
    summary={'cancer':'COAD','run_id':out.name,'version':spec['version'],'primary_n':33,'author_tumor_sensitivity_n':37,'families':{},'covariate_adjustment':'NOT_RUN','functional_evidence':'NOT_RUN'}
    for kind,rows in records.items():
        good=[r for r in rows if r['status']=='DONE'];qs=bh([r['p_value'] for r in good])
        for r,q in zip(good,qs):r['q_value']=float(q)
        for r in rows:
            if r['test_family']!='NOT_IN_CURRENT_FAMILY':r['family_n_evaluable']=len(good)
        rows.sort(key=lambda r:tuple(str(r[k] or '') for k in ['cohort','metabolite_name','metabolite_key','gene','human_gene_id']))
        assert len({(r['metabolite_name'],r['metabolite_key'],r['human_gene_id']) for r in rows})==len(rows)
        sig=[r for r in good if r['q_value']<.05]
        summary['families'][kind]={'planned':458 if kind=='paired_RNA33' else 674,'evaluable':len(good),'q_lt_005':len(sig),
          'significant_positive':sum(r['effect']>0 for r in sig),'significant_negative':sum(r['effect']<0 for r in sig),
          'significant_unique_genes':len({r['human_gene_id'] for r in sig}),
          'n_distribution':dict(Counter(str(r['n']) for r in good)),
          'not_evaluable_reasons':dict(Counter(r['reason'] for r in rows if r['status']=='NOT_EVALUABLE'))}
        table(out/('results.tsv' if kind=='primary33' else kind+'.tsv'),fields,rows)
    assert summary['families']['primary33']['evaluable']==spec['association_family_expected_evaluable']==652
    assert summary['families']['paired_RNA33']['evaluable']==spec['rna_family_expected_evaluable']==446
    assert all(r['n']==33 for r in records['primary33'] if r['status']=='DONE')
    # Every row carries exact frozen strings; unavailable results never acquire P/q.
    for kind in ['primary33','all37','availability33']:
        for r in records[kind]:
            f=effects[(r['metabolite_name'],r['metabolite_key'])]
            assert (r['original_effect'],r['original_q'])==(f['hedges_g'],f['effect_fdr'])
            if r['status']!='DONE':assert r['p_value'] is None and r['q_value'] is None
    for p,h in hashes.items():assert sha(p)==h,'Input changed during run'
    dump(out/'summary.json',summary)
    table(out/'source_manifest.tsv',['source_path','sha256'],[dict(source_path=str(p.relative_to(root)),sha256=h) for p,h in hashes.items()])
    outputs=['results.tsv','all37.tsv','availability33.tsv','paired_RNA33.tsv','summary.json','analysis_spec.json','source_manifest.tsv']
    dump(out/'validation.json',dict(status='PASS',checks=['locked family before calculation','explicit independent patients and exact paired controls','input hashes unchanged','all primary rho checked against scipy.spearmanr','all original effects/q preserved','family-specific BH','all unavailable P/q remain NA','unique relation keys','expected652 associations and446 paired RNA genes'],
      public_output_sha256={f:sha(out/f) for f in outputs},limitations=['No purity or batch adjustment','No covariate-adjusted inference in this batch','Author availability is not a verified detection mask','Within-CAMP exploratory analysis, not independent validation']))
    (out/'.running').unlink();(out/'DONE').write_text(datetime.now(timezone.utc).isoformat()+'\n')
    print(json.dumps(summary),flush=True)

if __name__=='__main__':
    try:main()
    except Exception:
        traceback.print_exc();sys.exit(1)
