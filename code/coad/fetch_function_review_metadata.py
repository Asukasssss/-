"""Fetch official bibliographic records and selected OA XML into a local cache."""
import argparse,concurrent.futures,hashlib,json,urllib.parse,urllib.request
from pathlib import Path
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--evidence-json',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=False);rows=json.loads(a.evidence_json.read_text(encoding='utf-8'))
    clauses=['(EXT_ID:'+r['PMID']+' AND SRC:MED)' if r['PMID'] else 'DOI:"'+r['DOI']+'"' for r in rows]
    url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':' OR '.join(clauses),'format':'json','resultType':'core','pageSize':100})
    with urllib.request.urlopen(url,timeout=50) as response:raw=response.read()
    (a.out/'europepmc_core.json').write_bytes(raw);(a.out/'request.json').write_text(json.dumps({'url':url,'sha256':hashlib.sha256(raw).hexdigest()},indent=2),encoding='utf-8')
    data=json.loads(raw)['resultList']['result']
    def fetch(ident):
        r=next(x for x in data if x['id']==ident);pmc=r['pmcid'];url='https://www.ebi.ac.uk/europepmc/webservices/rest/'+pmc+'/fullTextXML'
        try:
            with urllib.request.urlopen(url,timeout=45) as response:b=response.read()
            (a.out/(pmc+'.xml')).write_bytes(b);return dict(PMID=ident,pmcid=pmc,url=url,status='DOWNLOADED_FOR_SELECTED_PASSAGE_REVIEW',sha256=hashlib.sha256(b).hexdigest())
        except Exception as e:return dict(PMID=ident,pmcid=pmc,url=url,status='ACCESS_BLOCKED',error=type(e).__name__)
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:download=list(pool.map(fetch,['39819335','34539890','39642324','37343364','36688005']))
    (a.out/'fulltext_fetch.json').write_text(json.dumps(download,indent=2),encoding='utf-8')
    print(json.dumps({'metadata_records':len(data),'fulltext_fetches':download}))
if __name__=='__main__':main()
