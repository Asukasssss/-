"""Append source and paired-expression evidence to all39 without changing old fields."""
import argparse,csv,hashlib,json
from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parents[2]
RUN='20260920T091100Z_cell_paired35_v1'
SOURCE=R/'results/COAD/06_EXTERNAL/20260920T081200Z_cell_expression35_v2'

def read(p):
    with p.open(encoding='utf-8',newline='') as f:
        r=csv.DictReader(f,delimiter='\t');return r.fieldnames,list(r)

def direction(r):
    sign='T>N' if int(r['n_positive'])>int(r['n_negative']) else 'T<N'
    return sign+(';median_delta=0;nonzero='+r['n_informative']+'/'+r['n'] if float(r['effect'])==0 else '')

def main():
    d=R/'results/COAD/06_EXTERNAL'/RUN
    prior=R/'results/COAD/04_ROBUSTNESS/20260919T142500Z_covariates_v1/nominal39_with_covariates.tsv'
    header,old=read(prior);_,primary=read(d/'primary_results.tsv');_,source=read(SOURCE/'tumor_lineage_overview.tsv')
    gene_summary=[];full=[]
    genes=json.loads((d/'analysis_spec.json').read_text())['targets']
    for gene in genes:
        h,rows=read(d/'by_gene'/(gene+'.tsv'));full.extend(rows)
        assert len(rows)==766 and all(r['gene']==gene for r in rows)
        assert len({tuple(r[k] for k in ['study','enrichment','technology','annotation_level','cell_type']) for r in rows})==len(rows)
        for r in rows:
            if r['status']=='DONE':
                assert int(r['n'])>=6 and 0<=float(r['p_value'])<=float(r['q_value'])<=1
                assert float(r['ci_lower'])<=float(r['effect'])<=float(r['ci_upper']) and float(r['ci_coverage'])>=.95
                assert int(r['n_positive'])+int(r['n_negative'])+int(r['n_zero'])==int(r['n'])
            else:assert r['p_value']==r['q_value']=='NA'
            if int(r['n'])<3 or gene=='GSTT2':assert r['effect']=='NA'
        row=dict(gene=gene,full_strata=len(rows),missing_gene=gene=='GSTT2')
        for study in ['Lee','Pelka']:
            g=[r for r in rows if r['study'].startswith(study) and r['tier']=='primary'];test=[r for r in g if r['status']=='DONE'];sig=[r for r in test if float(r['q_value'])<.05]
            row.update({study+'_primary_evaluable':len(test),study+'_primary_p_lt05':sum(float(r['p_value'])<.05 for r in test),study+'_primary_q_lt05':len(sig),study+'_q_supported_contexts':'; '.join(r['cell_type']+':'+direction(r) for r in sig)or 'NONE'})
        gene_summary.append(row)
    assert len(full)==26810 and len(primary)==700
    # Independently recompute BH with scalar rank formulation, keeping all-zero tests.
    f=pd.DataFrame(full)
    for family,g in f[f.status=='DONE'].groupby('test_family'):
        ordered=sorted(g.to_dict('records'),key=lambda r:float(r['p_value']));m=len(ordered);running=1.
        for i in range(m-1,-1,-1):
            running=min(running,float(ordered[i]['p_value'])*m/(i+1));assert abs(running-float(ordered[i]['q_value']))<1e-12
            assert int(ordered[i]['family_n_evaluable'])==m
    pd.DataFrame(gene_summary).to_csv(d/'gene35_evidence_summary.tsv',sep='\t',index=False)
    new=[]
    fields=['study','enrichment','technology','cell_type','n','effect','ci_lower','ci_upper','p_value','q_value','status','reason','n_informative','paired_tumor_CPM_median','paired_normal_CPM_median']
    for original in old:
        gene=original['gene'];row=original.copy();g=[r for r in primary if r['gene']==gene]
        row.update(cell35_source_version='COAD_cell_origin35_v2',cell35_source_context_json=json.dumps([{k:r[k] for k in ['study','technology','highest_median_CPM_lineage','n_qualified_patients','median_CPM','median_detection']} for r in source if r['gene']==gene],ensure_ascii=False),
            cell35_paired_version='COAD_cell_paired35_v1',cell35_paired_primary_json=json.dumps([{k:r[k] for k in fields} for r in g],ensure_ascii=False),
            cell35_all_strata_path='results/COAD/06_EXTERNAL/'+RUN+'/by_gene/'+gene+'.tsv',
            cell35_metabolite_measurement='NOT_AVAILABLE',cell35_axis_validation='NOT_VALIDATED',
            cell35_decision='Preserve original9/21/5; external expression evidence is separate from CAMP metabolite-RNA association',
            cell35_remaining_gap='Independent matching metabolite-RNA validation; applicable perturbation/mediator evidence; source/measurement limits')
        assert all(row[k]==original[k] for k in header);new.append(row)
    assert len(new)==39 and len({r['gene'] for r in new})==35
    with (d/'relations39_integrated.tsv').open('w',encoding='utf-8',newline='') as out:
        w=csv.DictWriter(out,fieldnames=list(new[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(new)
    validation=json.loads((d/'validation.json').read_text());validation.update(local_validation='PASS',all39_original_fields_string_identical=True,
        independent_BH_recalculation=True,unique_keys=True,privacy_thresholds=True,missing_gene_retained=True,ci_bounds_and_coverage=True,
        source_relation_sha256=hashlib.sha256(prior.read_bytes()).hexdigest(),integration_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (d/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(pd.DataFrame(gene_summary).to_string(index=False))
if __name__=='__main__':main()
