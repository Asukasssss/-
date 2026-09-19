"""Validate public aggregate tables without accessing patient-level matrices."""
import argparse,csv,hashlib,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--run-id',required=True);a=p.parse_args()
    out=(ROOT/'results/COAD/03_PATIENT'/a.run_id).resolve();assert out.parent==(ROOT/'results/COAD/03_PATIENT').resolve()
    validation=json.loads((out/'validation.json').read_text());assert validation['status']=='PASS'
    for name,expected in validation['public_output_sha256'].items():assert digest(out/name)==expected
    spec=json.loads((out/'analysis_spec.json').read_text());assert spec['association_family_planned']==674
    template=(ROOT/'templates/statistical_result.tsv').read_text().strip().split('\t')
    effects={(r['feature_name'],r['metabolite_key']):r for r in read(ROOT/'reference/camp/cancer_effects.tsv') if r['cancer']=='COAD'}
    summary=json.loads((out/'summary.json').read_text())
    for kind in ['primary33','all37','availability33','paired_RNA33']:
        rows=read(out/('results.tsv' if kind=='primary33' else kind+'.tsv'))
        assert len(rows)==(458 if kind=='paired_RNA33' else 974)
        assert list(rows[0])[:len(template)]==template
        keys=[tuple(r[k] for k in ['cohort','metabolite_name','metabolite_key','gene','human_gene_id']) for r in rows]
        assert keys==sorted(keys) and len(keys)==len(set(keys))
        good=[r for r in rows if r['status']=='DONE']
        assert len(good)==summary['families'][kind]['evaluable']
        assert sum(float(r['q_value'])<.05 for r in good)==summary['families'][kind]['q_lt_005']
        ranked=sorted(good,key=lambda r:float(r['p_value']))
        tail=1.0
        for i in range(len(ranked)-1,-1,-1):
            tail=min(tail,len(ranked)*float(ranked[i]['p_value'])/(i+1))
            assert math.isclose(float(ranked[i]['q_value']),tail,rel_tol=1e-12,abs_tol=1e-14)
        for r in rows:
            if kind!='paired_RNA33':
                f=effects[(r['metabolite_name'],r['metabolite_key'])]
                assert (r['original_effect'],r['original_q'])==(f['hedges_g'],f['effect_fdr'])
            if r['status']=='DONE':
                assert 0<float(r['p_value'])<=float(r['q_value'])+1e-12<=1+1e-12
                assert int(r['family_n_evaluable'])==len(good)
                assert int(r['bootstrap_valid'])>=3800
                assert float(r['ci_lower'])<=float(r['ci_upper'])
                if kind!='paired_RNA33':
                    assert -1<=float(r['effect'])<=1
                    assert int(r['permutations'])==9999 and float(r['p_value'])>=.0001
                    assert -1.00000001<=float(r['ci_lower'])<=float(r['ci_upper'])<=1.00000001
                if kind in ['primary33','paired_RNA33']:assert r['n']=='33'
                if kind=='all37':assert r['n']=='37'
                if kind=='availability33':assert 8<=int(r['n'])<=33
                if kind=='primary33':assert r['loo_valid']=='33'
            else:
                assert r['p_value']==r['q_value']==r['effect']=='NA'
    headers=[x[3:] for x in (out/'README_CN.md').read_text(encoding='utf-8').splitlines() if x.startswith('## ')]
    assert headers==['本轮问题','输入与范围','实际结果','新手解释','限制/反证','当前决定','下一步','复现命令']
    index=read(ROOT/'coordination/stages/COAD.tsv')
    assert [(r['stage_id'],r['run_id']) for r in index]==sorted((r['stage_id'],r['run_id']) for r in index)
    print('PASS: public hashes, common fields/order/unique keys, frozen values, sample sizes, resampling coverage, all family-specific BH values and missing-result states.')

if __name__=='__main__':main()
