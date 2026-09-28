"""Bounded alias follow-up for zero-hit candidates and PNP name ambiguity."""
import csv,json,hashlib,time,urllib.request,urllib.parse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/PDAC/05_FUNCTION/20260920T114000Z_targeted_followup_v1'
CACHE=ROOT/'runtime/pdac_function_alias_v1';CACHE.mkdir(parents=True,exist_ok=True)
def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter='\t'))
def write(p,rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def main():
    initial=read(ROOT/'results/PDAC/05_FUNCTION/20260920T112000Z_all_gene_search_v1/gene_search_ledger.tsv')
    aliases=json.loads((ROOT/'results/PDAC/02_MAPPING/20260920T111500Z_direct_mapping_v1/gene_aliases.json').read_text())
    targets={r['gene']:[x for x in aliases[r['gene']] if len(x)>=3][:3] for r in initial if r['hit_count']=='0'}
    targets['PNP']=['purine nucleoside phosphorylase']
    def get(item):
        gene,names=item
        row={'gene':gene,'aliases':';'.join(names) or 'NA','query':'NA','status':'NOT_EVALUABLE','hits':'NA','response_sha256':'NA','retrieved_utc':datetime.now(timezone.utc).isoformat(),'reason':'No reviewed aliases available;absence of hits is not no function'}
        papers=[]
        if not names:return row,papers
        query='(TITLE_ABS:pancreatic OR TITLE_ABS:PDAC) AND ('+' OR '.join('TITLE_ABS:"'+n+'"' for n in names)+')'
        url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':query,'format':'json','resultType':'core','pageSize':3})
        row['query']=url;p=CACHE/(gene+'.json')
        try:
            if p.exists():data=p.read_bytes()
            else:
                time.sleep(.4)
                with urllib.request.urlopen(url,timeout=25) as f:data=f.read()
                p.write_bytes(data)
            obj=json.loads(data);row.update(status='DONE',hits=obj['hitCount'],response_sha256=hashlib.sha256(data).hexdigest(),reason='Alias retrieval only;not verified gene-specific intervention')
            for r in obj.get('resultList',{}).get('result',[]):papers.append({'gene':gene,'pmid':r.get('pmid','NA'),'article_id':r['id'],'source':r['source'],'title':r.get('title','NA'),'url':'https://europepmc.org/article/'+r['source']+'/'+r['id'],'status':'LOCATED_NOT_ADJUDICATED'})
        except Exception as e:row.update(status='ACCESS_BLOCKED',reason=str(e))
        return row,papers
    with ThreadPoolExecutor(max_workers=3) as pool:results=list(pool.map(get,sorted(targets.items())))
    rows=[r for r,p in results];papers=[p for r,ps in results for p in ps]
    write(OUT/'alias_search_ledger.tsv',rows)
    if papers:write(OUT/'alias_located_articles.tsv',papers)
    summary={'zero_hit_genes_followed':len(targets)-1,'queries_with_aliases':sum(r['status']=='DONE' for r in rows),'followup_genes_with_hits':sum(isinstance(r['hits'],int) and r['hits']>0 for r in rows),'located_records':len(papers),'PNP_full_name_followup':next(r for r in rows if r['gene']=='PNP')}
    (OUT/'alias_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary))
if __name__=='__main__':main()
