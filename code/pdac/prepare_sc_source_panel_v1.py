"""Preserve all51 paired nominal hits and existing mapping limits; no new biochemistry."""
import csv,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
RUN='20260921T151600Z_sc_source_three_cohorts_v1'
OUT=ROOT/'results/PDAC/06_EXTERNAL'/RUN
def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
def write(p,rows):
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def main():
 OUT.mkdir(parents=True,exist_ok=False)
 paired=ROOT/'results/PDAC/04_ROBUSTNESS/20260921T121509Z_paired_metabolites_v1/results.tsv'
 mapping=ROOT/'results/PDAC/02_MAPPING/20260920T111500Z_direct_mapping_v1'
 hits=[r for r in read(paired) if r['analysis_type']=='primary' and float(r['p_value'])<.05]
 assert len(hits)==51
 direct=read(mapping/'direct_relations_v1.tsv');conditional=read(mapping/'conditional_relations_v1.tsv')
 ledger=[];genes=set();relations=[]
 for r in hits:
  name=r['metabolite_name'];a=[e for e in direct if e['feature_name']==name];b=[e for e in conditional if e['feature_name']==name]
  gs=sorted({e['gene'] for e in a+b});genes.update(gs)
  ledger.append({'feature_name':name,'paired_p':r['p_value'],'paired_q':r['q_value'],'direct_gene_n':len({e['gene'] for e in a}),'conditional_gene_n':len({e['gene'] for e in b}),'genes':';'.join(gs),'mapping_status':'EXISTING_DIRECT_AND_OR_CONDITIONAL' if gs else 'MAPPING_NOT_YET_ESTABLISHED','mapping_limit':'Existing old42-derived ledger only;no new identity or exact-substrate certification'})
  for group,kind in [(a,'DIRECT_ANNOTATION'),(b,'CONDITIONAL')]:
   for e in group:relations.append({'feature_name':name,'gene':e['gene'],'mapping_class':kind,'source_url':e.get('source_urls',e.get('source_url','')),'limitation':e.get('reason',e.get('identity_note',''))})
 write(OUT/'paired51_mapping_coverage.tsv',ledger);write(OUT/'relation_scope.tsv',relations)
 (OUT/'panel.json').write_text(json.dumps({'genes':sorted(genes),'scope':'Existing direct/conditional mapping of paired51 nominal hits;unmapped features retained','n_features':51,'mapped_features':sum(bool(r['genes']) for r in ledger),'unmapped_features':[r['feature_name'] for r in ledger if not r['genes']]},indent=2)+'\n')
 sources=[paired,mapping/'direct_relations_v1.tsv',mapping/'conditional_relations_v1.tsv']
 write(OUT/'panel_manifest.tsv',[{'path':str(p.relative_to(ROOT)).replace('\\','/'),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources])
 print(json.dumps({'genes':len(genes),'features_with_mapping':sum(bool(r['genes']) for r in ledger),'unmapped':sum(not r['genes'] for r in ledger)}))
if __name__=='__main__':main()
