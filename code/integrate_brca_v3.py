"""Aggregate-only evidence update. Old functional columns are retained explicitly."""
import argparse,json
from pathlib import Path
import pandas as pd

ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);a=ap.parse_args();p=a.results
d=pd.read_csv(p/'candidate_comparison_117_v3.tsv',sep='\t');r=pd.read_csv(p/'er_availability_all_174.tsv',sep='\t');studies=json.loads((p/'functional_review_update.json').read_text())
assert len(d)==117 and d.gene.is_unique and len(r)==174 and r.relation_id.is_unique
both=r[r.primary_er_q.lt(.05)&r.q_value.lt(.05)]
d['n_same_relations_supported_both_ER_families']=d.gene.map(both.groupby('gene').size()).fillna(0).astype(int)
d['functional_review_this_round']='NOT_REVIEWED_THIS_ROUND; previous evidence retained'
d['functional_update_sources']='';d['functional_update_conclusion']='';d['next_action']='Retain pool; await suitable dependency data or manual functional review'
for gene in {s['gene'] for s in studies}:
    selected=[s for s in studies if s['gene']==gene];ix=d.gene.eq(gene)
    d.loc[ix,'functional_review_this_round']='; '.join(sorted({s['review_scope'] for s in selected}))
    d.loc[ix,'functional_update_sources']='; '.join(s['source_url'] for s in selected)
    d.loc[ix,'functional_update_conclusion']=' | '.join(str(s['year'])+': '+s['decision']+' '+s['limits'] for s in selected)
    d.loc[ix,'next_action']='Prioritize context-specific followup; evaluate research increment and metabolite mediation separately'
d.loc[d.gene.eq('GPCPD1'),'novelty_notes']='Existing breast cancer intervention and GPC mechanistic work in 2023/2024; CAMP correlation adds patient context, not first target discovery.'
d['depmap_status']='ACCESS_BLOCKED_HTTP403_RECHECKED_20260919'
d['followup_status']='Availability-subset ER sensitivity completed; functional review updated for GPCPD1/GPI only; all 117 genes retained'
d.to_csv(p/'candidate_comparison_117_integrated_v3.tsv',sep='\t',index=False)
pd.DataFrame(studies).to_csv(p/'functional_review_update.tsv',sep='\t',index=False)
print('Integrated',len(d),'genes; genes supported by same relationships in both ER families:',d.loc[d.n_same_relations_supported_both_ER_families.gt(0),'gene'].tolist())
