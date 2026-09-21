"""Bounded public literature discovery. Hits are not functional validation."""
from pathlib import Path
import sys,json,re,hashlib,concurrent.futures,time
import requests,pandas as pd
R=Path(sys.argv[1]);repo=Path(__file__).resolve().parents[1];R.mkdir(parents=True,exist_ok=True);(R/'sources').mkdir(exist_ok=True)
genes=pd.read_csv(repo/'results/BRCA/02_MAPPING/20260921T124136Z_mapping190_v2/genes_unique.tsv',sep='\t');genes=genes.loc[~genes.in_original117,'gene'].tolist();assert len(genes)==39
focused='--focused' in sys.argv
aliases={'MGLL':'(MGLL OR MAGL OR "monoacylglycerol lipase")','ENPP2':'(ENPP2 OR autotaxin)','PTGS1':'(PTGS1 OR COX1 OR "COX-1")','PTGS2':'(PTGS2 OR COX2 OR "COX-2")','LYPLA1':'(LYPLA1 OR APT1)','LYPLA2':'(LYPLA2 OR APT2)','LPCAT4':'(LPCAT4 OR AGPAT7 OR LPEAT2)','PUDP':'(PUDP OR HDHD1)','SLC22A2':'(SLC22A2 OR OCT2)','SLC22A4':'(SLC22A4 OR OCTN1)','SLC47A1':'(SLC47A1 OR MATE1)','SLC47A2':'(SLC47A2 OR MATE2)','ALOX15':'(ALOX15 OR "15-lipoxygenase")','ALOX5':'(ALOX5 OR "5-lipoxygenase")'}
def one(g):
 term='(LPCAT4 OR AGPAT7 OR LPEAT2)' if g=='LPCAT4' else ('(PUDP OR HDHD1)' if g=='PUDP' else g)
 query=f'TITLE_ABS:{term} AND (TITLE_ABS:"breast cancer" OR TITLE_ABS:"mammary carcinoma" OR TITLE_ABS:"breast carcinoma") AND (knockdown OR knockout OR silencing OR depletion OR overexpression OR CRISPR OR inhibition)'
 if focused:query=f'TITLE:{aliases.get(g,g)} AND ("breast cancer" OR "mammary" OR "breast carcinoma") NOT PUB_TYPE:Review'
 p=R/'sources'/(g+('_focused' if focused else '')+'.json');row=dict(gene=g,query=query,endpoint='https://www.ebi.ac.uk/europepmc/webservices/rest/search',retrieved_at_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),status='DONE')
 try:
  if not p.exists() or 'hitCount' not in json.loads(p.read_text(encoding='utf-8')):
   r=requests.get(row['endpoint'],params={'query':query,'format':'json','resultType':'core','pageSize':20},timeout=45);r.raise_for_status();p.write_text(json.dumps(r.json(),ensure_ascii=False),encoding='utf-8')
  j=json.loads(p.read_text(encoding='utf-8'));row.update(hit_count=j['hitCount'],returned=len(j['resultList']['result']),sha256=hashlib.sha256(p.read_bytes()).hexdigest());out=[]
  for x in j['resultList']['result']:
   abstract=re.sub('<[^>]+>',' ',x.get('abstractText',''));title=x.get('title','');url=('https://pubmed.ncbi.nlm.nih.gov/'+x['pmid']+'/') if x.get('pmid') else ('https://doi.org/'+x['doi'] if x.get('doi') else '')
   out.append(dict(gene=g,pmid=x.get('pmid',''),pmcid=x.get('pmcid',''),doi=x.get('doi',''),title=title,year=x.get('pubYear',''),abstract=abstract,url=url,publication_types=';'.join(x.get('pubTypeList',{}).get('pubType',[])),gene_in_title=bool(re.search(r'\b'+g+r'\b',title,re.I)),retrieval_depth='title_and_abstract',functional_status='NEEDS_REVIEW',reason='search_hit_not_functional_evidence',geo_accessions_in_abstract=';'.join(sorted(set(re.findall(r'GSE\d+',abstract))))))
  return row,out
 except Exception as e:row.update(status='ACCESS_BLOCKED',reason=type(e).__name__+':'+str(e)[:150]);return row,[]
rows=[];hits=[]
for row,out in concurrent.futures.ThreadPoolExecutor(max_workers=3).map(one,genes):rows.append(row);hits+=out;print(row['gene'],row.get('hit_count',row['status']),flush=True)
suffix='_focused' if focused else ''
pd.DataFrame(rows).to_csv(R/('search_log'+suffix+'.tsv'),sep='\t',index=False);pd.DataFrame(hits).to_csv(R/('literature_hits'+suffix+'.tsv'),sep='\t',index=False)
(R/'analysis_spec.json').write_text(json.dumps(dict(version='new39_literature_v1',genes=genes,source_commit='9e4387f',query_scope='gene plus breast cancer plus intervention vocabulary;top20 relevance hits;title+abstract only',not_systematic_review=True,hits_not_functional_validation=True,LPCAT4_requires_disambiguation='Q643R3 LPCAT4/AGPAT7/LPEAT2 vs Q6ZWT7 MBOAT2 historical LPCAT4',no_DepMap=True),indent=2))
