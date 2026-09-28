"""Targeted patient-paired comparisons; cell/patient values stay on server165."""
import argparse, hashlib, json, pathlib, platform
import numpy as np
import pandas as pd
import scipy
from scipy.stats import wilcoxon

COMPS=['Myeloid','Lymphoid','Oligodendrocyte','Astrocyte','Neuron','OPC','Vascular','All_nonmalignant']
PREFIX='cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA')
def bh(p):
    p=np.asarray(p);o=np.argsort(p);v=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];a=np.empty(len(p));a[o]=np.minimum(v,1);return a
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def main(root,old):
    out=root/'public';rng=np.random.default_rng(20260928)
    snu=pd.read_csv(root/'private/SNUH_patient_profiles.tsv',sep='\t')
    care=pd.read_csv(old/'private/sc_donor_profiles_private.tsv',sep='\t')
    care=care[care.gene.eq('LYPLA1')].rename(columns={'Patient ID':'patient','celltype':'category'})
    care['category']=care.category.replace({'Neoplastic':'Malignant'});care['cohort']='CARE10x'
    cols=['cohort','partition','patient','category','n_cells','sum_expression','sum_detected']
    d=pd.concat([snu[cols],care[cols]],ignore_index=True)
    assert not d.duplicated(['cohort','partition','patient','category']).any()
    nonmal=d[~d.category.isin(['Malignant','Unresolved','Other'])]
    pool=nonmal.groupby(['cohort','partition','patient'])[['n_cells','sum_expression','sum_detected']].sum().reset_index()
    pool['category']='All_nonmalignant';d=pd.concat([d,pool],ignore_index=True)
    d['mean_expression']=d.sum_expression/d.n_cells;d['detection']=d.sum_detected/d.n_cells
    save(d,root/'private/combined_patient_profiles.tsv')
    result=[];profiles=[]
    for (cohort,partition),part in d.groupby(['cohort','partition'],sort=True):
        scale='author log2(CPM+1)' if cohort=='SNUH' else 'log1p(CP10k)'
        for threshold in [20,50]:
            z=part[part.n_cells>=threshold]
            for category,t in z.groupby('category'):
                profiles.append(dict(cohort=cohort,partition=partition,min_cells=threshold,category=category,
                    patients=len(t),cells=int(t.n_cells.sum()),mean_expression=t.mean_expression.mean(),
                    median_expression=t.mean_expression.median(),mean_detection=t.detection.mean(),scale=scale))
            family=f'{cohort}:{partition}:min{threshold}:LYPLA1:8_predefined_comparators'
            rows=[]
            for comparator in COMPS:
                a=z[z.category.eq('Malignant')].set_index('patient');b=z[z.category.eq(comparator)].set_index('patient')
                ids=a.index.intersection(b.index);av=a.loc[ids,'mean_expression'];bv=b.loc[ids,'mean_expression'];delta=(av-bv).to_numpy()
                n=len(ids);good=n>=10
                row=dict(cancer='GBM',cohort=cohort,stage_id='06_EXTERNAL',run_id=root.name,analysis_version='gbm_lypla1_malignant_v3',
                    analysis_type='paired_patient_expression',metabolite_key='NA',metabolite_name='NA',gene='LYPLA1',unit='patient',n=n,n_reference=n,
                    effect_type='mean paired difference:Malignant-minus-comparator',effect=float(delta.mean()) if n else np.nan,
                    ci_lower=np.nan,ci_upper=np.nan,p_value=np.nan,q_value=np.nan,test_family=family,family_n_evaluable=0,
                    status='DONE' if good else 'NOT_EVALUABLE',reason='targeted_postselection_inference' if good else 'fewer_than_10_matched_patients',
                    source_id='SNUH_2026_UCSC' if cohort=='SNUH' else 'CARE_2025_GSE274546',partition=partition,comparator=comparator,
                    min_cells=threshold,effect_scale=scale,positive_patients=int((delta>0).sum()) if n else np.nan,
                    positive_fraction=float((delta>0).mean()) if n else np.nan,median_paired_difference=float(np.median(delta)) if n else np.nan,
                    malignant_mean=float(av.mean()) if n else np.nan,comparator_mean=float(bv.mean()) if n else np.nan,
                    malignant_detection=float(a.loc[ids,'detection'].mean()) if n else np.nan,
                    comparator_detection=float(b.loc[ids,'detection'].mean()) if n else np.nan,
                    malignant_cells=int(a.loc[ids,'n_cells'].sum()) if n else 0,comparator_cells=int(b.loc[ids,'n_cells'].sum()) if n else 0,
                    family_n_predefined=8,bootstrap_draws=10000 if good else 0)
                if good:
                    row['p_value']=1.0 if np.all(delta==0) else float(wilcoxon(delta,alternative='two-sided',zero_method='wilcox',method='auto').pvalue)
                    draws=delta[rng.integers(0,n,size=(10000,n))].mean(axis=1)
                    row['ci_lower'],row['ci_upper']=np.quantile(draws,[.025,.975])
                rows.append(row)
            inds=[i for i,r in enumerate(rows) if r['status']=='DONE']
            q=bh([rows[i]['p_value'] for i in inds])
            for r in rows:r['family_n_evaluable']=len(inds)
            for i,v in zip(inds,q):rows[i]['q_value']=v
            result.extend(rows)
    r=pd.DataFrame(result);r=r[PREFIX+[c for c in r if c not in PREFIX]];save(r,out/'LYPLA1_paired_results.tsv')
    save(pd.DataFrame(profiles),out/'LYPLA1_celltype_profiles.tsv')
    spec=json.loads((root/'source/analysis_spec_pre_expression.json').read_text());spec['cohort_selection']=json.loads((root/'source/cohort_frozen_before_expression.json').read_text())
    spec['software']=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__)
    spec['normalization']=json.loads((root/'source/normalization_verified.json').read_text())
    spec['CARE_all_nonmalignant']='known nonmalignant categories only;Other/Unresolved excluded;pool all cells before threshold'
    (out/'analysis_spec.json').write_text(json.dumps(spec,indent=2))
    paths=[p for p in (root/'source').iterdir() if p.is_file() and not p.name.startswith('meta.part.') and p.name not in ['SNUH_meta.tsv','SNUH_meta_parallel.tsv']]
    paths += list((root/'code').glob('*'))+[old/'private/sc_donor_profiles_private.tsv']
    save(pd.DataFrame([dict(source_path=str(p),bytes=p.stat().st_size,sha256=sha(p),scope='server_only_input_or_code') for p in paths]),out/'source_manifest.tsv')
    main=r[r.min_cells.eq(20)&((r.cohort.eq('SNUH')&r.partition.eq('all'))|r.cohort.eq('CARE10x'))]
    print(main[['cohort','partition','comparator','n','effect','ci_lower','ci_upper','p_value','q_value','positive_fraction','status']].to_string(index=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=pathlib.Path);p.add_argument('old',type=pathlib.Path);a=p.parse_args();main(a.root,a.old)
