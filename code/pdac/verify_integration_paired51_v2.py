from pathlib import Path
import json
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
def main():
 p=ROOT/'results/PDAC/07_INTEGRATION/20260922T055500Z_internal_sc_paired51_v2';a=pd.read_csv(p/'candidate_pre_scRNA.tsv',sep='\t').set_index('relation_id');b=pd.read_csv(p/'candidate_scRNA_appended.tsv',sep='\t').set_index('relation_id');pd.testing.assert_frame_equal(a,b[a.columns],check_exact=True)
 source=pd.read_csv(ROOT/'results/PDAC/06_EXTERNAL/20260922T054000Z_sc_paired51_v2/cross_cohort_source.tsv',sep='\t').set_index('gene')
 allfields=source.columns.tolist();actual=b[b.gene!='UNRESOLVED'].set_index('gene')[allfields];expected=source.loc[actual.index];pd.testing.assert_frame_equal(actual,expected,check_exact=False,check_dtype=False,rtol=0,atol=1e-12)
 numeric=expected.select_dtypes(include='number').columns;max_delta=float((actual[numeric]-expected[numeric]).abs().max().max())
 assert len(b)==722 and b.metabolite_name.nunique()==51 and len(set(b.gene)-{'UNRESOLVED'})==538 and len(b[b.mapping_tier!='UNRESOLVED'])==715
 assert b.metabolite_paired_p.lt(.05).all() and not b.metabolite_paired_q.lt(.05).any()
 assert b.exact_relation_external_validation.eq('NOT_PERFORMED').all() and b.functional_perturbation.eq('NOT_PERFORMED').all()
 summary=json.loads((p/'summary.json').read_text());assert summary['candidate_classes']=={x:int(b.candidate_class.eq(x).sum()) for x in 'ABCD'}
 result={'status':'PASS','all722_internal_rows_and_columns_unchanged_after_scRNA_append':True,'all715_relation_source_joins_match_gene_source_table':True,'source_numeric_csv_roundtrip_tolerance':1e-12,'max_observed_source_numeric_abs_delta':max_delta,'unresolved_features_retained':7,'source_gene_count':538,'all51_nominal_not_FDR_supported':True,'candidate_class_counts_reconciled':True}
 (p/'independent_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
