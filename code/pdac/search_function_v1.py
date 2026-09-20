"""One bounded, logged literature search for every locked PDAC candidate.
Text flags are retrieval triage, never confirmation of a genetic intervention.
"""
import csv,hashlib,json,re,time,urllib.request,urllib.parse
from pathlib import Path
from datetime import datetime,timezone
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[2]
MAP=ROOT/'results/PDAC/02_MAPPING/20260920T111500Z_direct_mapping_v1'
CACHE=ROOT/'runtime/pdac_function_v1';CACHE.mkdir(parents=True,exist_ok=True)
OUT=ROOT/'results/PDAC/05_FUNCTION/20260920T112000Z_all_gene_search_v1'
def write(p,rows,fields=None):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def get(g):
    query=f'(TITLE_ABS:pancreatic OR TITLE_ABS:PDAC) AND TITLE_ABS:"{g}"'
    url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':query,'format':'json','resultType':'core','pageSize':5})
    p=CACHE/(g+'.json');err='';data=None
    for retry in range(2):
        try:
            if p.exists():raw=p.read_bytes()
            else:
                time.sleep(.4)
                with urllib.request.urlopen(url,timeout=25) as f:raw=f.read()
                p.write_bytes(raw)
            data=json.loads(raw);break
        except Exception as e:err=type(e).__name__+': '+str(e);time.sleep(1)
    row={'gene':g,'query':query,'search_url':url,'search_status':'DONE' if data is not None else 'ACCESS_BLOCKED',
      'hit_count':data.get('hitCount',0) if data is not None else 'NA','records_located':len(data.get('resultList',{}).get('result',[])) if data is not None else 0,
      'review_status':'TITLE_ABSTRACT_RETRIEVAL_ONLY_NOT_FUNCTIONAL_CONFIRMATION','pdac_genetic_intervention_verified':'NOT_REVIEWED',
      'exact_metabolite_mechanism_verified':'NOT_REVIEWED','model_applicability':'NOT_REVIEWED','negative_or_opposite_evidence':'NOT_REVIEWED',
      'read_depth':'Metadata and abstract retrieved; keyword flags only','retrieved_utc':datetime.now(timezone.utc).isoformat(),
      'response_sha256':hashlib.sha256(raw).hexdigest() if data is not None else 'NA', 'error':err if data is None else 'NA'}
    papers=[]
    if data is not None:
        for d in data.get('resultList',{}).get('result',[]):
            text=d.get('title','')+' '+re.sub('<[^>]+>',' ',d.get('abstractText',''))
            papers.append({'gene':g,'source':d.get('source',''),'article_id':d.get('id',''),'pmid':d.get('pmid',''),'pmcid':d.get('pmcid',''),
              'doi':d.get('doi',''),'title':d.get('title',''),'year':d.get('pubYear',''),
              'url':'https://europepmc.org/article/'+d.get('source','MED')+'/'+d.get('id',''),
              'text_mentions_perturbation':bool(re.search(r'knockdown|knockout|silenc|overexpress|CRISPR|delet|shRNA|siRNA',text,re.I)),
              'text_mentions_metabolism':bool(re.search(r'metaboli|flux|isotope|rescue',text,re.I)),
              'text_mentions_PDCA_model':bool(re.search(r'adenocarcinoma|PDAC|PANC-1|MIA PaCa|BxPC|organoid|KPC',text,re.I)),
              'evidence_status':'LOCATED_NOT_VALIDATED','intervention_direction':'NOT_REVIEWED','metabolite_or_flux_measured':'NOT_REVIEWED',
              'rescue':'NOT_REVIEWED','opposing_result':'NOT_REVIEWED','reading_depth':'ABSTRACT_RETRIEVED_NOT_MANUALLY_ADJUDICATED'})
    return row,papers
def main():
    relations=list(csv.DictReader((MAP/'direct_relations_v1.tsv').open(encoding='utf-8'),delimiter='\t'))
    genes=sorted({r['gene'] for r in relations});rows=[];papers=[]
    OUT.mkdir(parents=True,exist_ok=True)
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(get,g) for g in genes]
        for i,f in enumerate(as_completed(futures),1):
            r,pp=f.result();rows.append(r);papers.extend(pp)
            if i%20==0:print('searched',i,'/',len(genes),flush=True)
    rows.sort(key=lambda r:r['gene']);papers.sort(key=lambda r:(r['gene'],r['article_id']))
    write(OUT/'gene_search_ledger.tsv',rows)
    if papers:write(OUT/'located_articles.tsv',papers)
    summary={'planned_genes':len(genes),'searched_successfully':sum(r['search_status']=='DONE' for r in rows),
      'blocked':sum(r['search_status']=='ACCESS_BLOCKED' for r in rows),'genes_with_hits':sum(isinstance(r['hit_count'],int) and r['hit_count']>0 for r in rows),
      'located_records':len(papers),'functional_review_status':'PARTIAL_RETRIEVAL_ONLY','verified_functional_genes':'NOT_YET_ADJUDICATED'}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
