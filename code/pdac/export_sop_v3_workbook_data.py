"""Human-facing workbook inputs; retain precision, all relations and all genes."""
import json,math
import pandas as pd
from record_sop_v3_resume import REPO,INTERNAL,DISCOVERY,INTEGRATION

def records(df):
    return [{k:(None if isinstance(v,float) and not math.isfinite(v) else v) for k,v in r.items()} for r in df.to_dict('records')]
def main():
    out=REPO/'outputs/pdac-sop-v3-20260922';out.mkdir(parents=True,exist_ok=True)
    rel=pd.read_csv(INTEGRATION/'candidate_relations_integrated.tsv',sep='\t');rel['order']=rel.mapping_status.map({'DIRECT':0,'CONDITIONAL':1});rel=rel.sort_values(['order','association_q','metabolite_name','gene']).drop(columns='order')
    genes=pd.read_csv(INTEGRATION/'candidate_genes_integrated.tsv',sep='\t').sort_values('gene');rna=pd.read_csv(INTERNAL/'RNA_P005_view.tsv',sep='\t').sort_values(['q_value','gene']);met=pd.read_csv(DISCOVERY/'metabolite_workpool.tsv',sep='\t').sort_values(['p_value','metabolite_name']);av=pd.read_csv(DISCOVERY/'metabolite_available_sensitivity.tsv',sep='\t').set_index('metabolite_key')
    for col in ['n','p_value','q_value','effect']:met['available_'+col]=av.loc[met.metabolite_key,col].to_numpy()
    assert len(rel)==715 and len(genes)==687 and len(rna)==91 and len(met)==51
    d={'relations':records(rel),'genes':records(genes),'RNA':records(rna),'metabolites':records(met),'summary':json.loads((INTEGRATION/'summary.json').read_text())}
    (out/'workbook_data.json').write_text(json.dumps(d,ensure_ascii=False,allow_nan=False),encoding='utf8');print(out)
if __name__=='__main__':main()
