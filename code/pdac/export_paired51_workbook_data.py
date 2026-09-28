"""Prepare aggregate-only scientific workbook views; underlying full TSVs remain authoritative."""
from pathlib import Path
import json
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
def clean(df):return df.astype(object).where(pd.notna(df),None).to_dict(orient='records')
def main():
 out=ROOT/'outputs/pdac-internal-sc-20260922';out.mkdir(parents=True,exist_ok=True)
 base=ROOT/'results/PDAC/07_INTEGRATION/20260922T055500Z_internal_sc_paired51_v2'
 c=pd.read_csv(base/'candidate_scRNA_appended.tsv',sep='\t');c['_tier_order']=c.mapping_tier.map({'DIRECT':0,'CONDITIONAL':1,'UNRESOLVED':2});c=c.sort_values(['_tier_order','CAMP_q','CAMP_p','metabolite_name','gene'],na_position='last').drop(columns='_tier_order')
 features=c.drop_duplicates('metabolite_name').sort_values('metabolite_name')
 genes=pd.read_csv(ROOT/'results/PDAC/03_PATIENT/20260922T053500Z_internal_paired51_v2/paired_RNA.tsv',sep='\t').merge(pd.read_csv(ROOT/'results/PDAC/06_EXTERNAL/20260922T054000Z_sc_paired51_v2/cross_cohort_source.tsv',sep='\t'),on='gene',validate='one_to_one')
 summary=json.loads((base/'summary.json').read_text());(out/'workbook_data.json').write_text(json.dumps({'relations':clean(c),'features':clean(features),'genes':clean(genes),'summary':summary},ensure_ascii=False,indent=2),encoding='utf-8')
if __name__=='__main__':main()
