"""Export original COAD nominal P view without recomputing any statistic."""
import argparse,csv,hashlib,json
from pathlib import Path


def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--frozen',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--template',type=Path,required=True);a=ap.parse_args()
    assert a.out.is_dir() and not (a.out/'results.tsv').exists()
    frozen={r['feature_name']:r for r in read(a.frozen) if r['cancer']=='COAD'}
    rows=[r for r in read(a.source) if r['cancer']=='COAD'];assert len(rows)==len(frozen)==159
    assert len({r['feature_name'] for r in rows})==159
    for r in rows:
        f=frozen[r['feature_name']]
        assert r['hedges_g']==f['hedges_g'] and r['wilcoxon_fdr']==f['effect_fdr']
        assert r['metabolite_key']==f['metabolite_key']
        assert r['n_tumor']=='37' and r['n_normal']=='39' and r['analysis_design'].startswith('Unpaired;')
    with (a.out/'original_coad159.tsv').open('x',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    selected=[r for r in rows if float(r['wilcoxon_p'])<.05]
    prefix=a.template.read_text().strip().split('\t')
    fields=prefix+['original_P_lt_005','original_q_lt_005','original_analysis_design']
    with (a.out/'results.tsv').open('x',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for r in rows:
            x=dict.fromkeys(fields,'NA')
            x.update(cancer='COAD',cohort='COAD',stage_id='01_CAMP',run_id=a.out.name,analysis_version='original_camp_p_v1',
                analysis_type='original_unpaired_P_view_no_recalculation',metabolite_key=r['metabolite_key'],metabolite_name=r['feature_name'],
                unit='specimen_original_unpaired_groups',n=r['n_tumor'],n_reference=r['n_normal'],effect_type='Hedges_g',effect=r['hedges_g'],
                ci_lower=r['hedges_g_ci_lower'],ci_upper=r['hedges_g_ci_upper'],p_value=r['wilcoxon_p'],q_value=r['wilcoxon_fdr'],
                test_family='Original COAD Wilcoxon family; preserved without recalculation',family_n_evaluable=159,status='DONE',
                reason='ORIGINAL_STATISTICS_PRESERVED',source_id='cohort_effects.tsv',original_P_lt_005=str(float(r['wilcoxon_p'])<.05),
                original_q_lt_005=str(float(r['wilcoxon_fdr'])<.05),original_analysis_design=r['analysis_design'])
            w.writerow(x)
    summary=dict(total=159,p_lt_005=len(selected),q_lt_005=sum(float(r['wilcoxon_fdr'])<.05 for r in rows),
        p_up=sum(float(r['hedges_g'])>0 for r in selected),p_down=sum(float(r['hedges_g'])<0 for r in selected),
        nominal_only=sum(float(r['wilcoxon_fdr'])>=.05 for r in selected),n_tumor=37,n_normal=39,new_tests=False)
    for name,obj in [('summary.json',summary),('analysis_spec.json',dict(version='original_camp_p_v1',selection='original wilcoxon_p < 0.05',scope='all159 frozen COAD',recalculate=False,grouping='sign of original Hedges g',unit='original unpaired specimens;37 tumor,39 normal',frozen_validation='All159 original g/q/metabolite keys exact strings match'))]:
        (a.out/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    with (a.out/'source_manifest.tsv').open('x',encoding='utf-8',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['source_path','version','sha256'])
        for p in [a.source,a.frozen,Path(__file__).resolve()]:w.writerow([str(p),'read_only_at_export',sha(p)])
    (a.out/'validation.json').write_text(json.dumps(dict(status='PASS',all159_original_g_q_keys_exact=True,original_rows_preserved=True,
        new_tests=False,output_sha256={n:sha(a.out/n) for n in ['original_coad159.tsv','results.tsv']}),indent=2)+'\n')
    print(json.dumps(summary))


if __name__=='__main__':main()
