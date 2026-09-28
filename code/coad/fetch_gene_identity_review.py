"""Fetch authoritative identity records for missing COAD RNA symbols; no matrix input."""
import csv, hashlib, json, time, urllib.parse, urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
CACHE=ROOT/'runtime/coad_gene_identity_20260919'


def main():
    CACHE.mkdir(parents=True,exist_ok=True)
    reg_path=CACHE/'source_registry.json'
    registry=json.loads(reg_path.read_text()) if reg_path.exists() else {}
    def fetch(url):
        path=CACHE/(hashlib.sha256(url.encode()).hexdigest()+'.json')
        if url in registry and path.exists():
            data=path.read_bytes();assert hashlib.sha256(data).hexdigest()==registry[url]['sha256']
            return json.loads(data)
        request=urllib.request.Request(url,headers={'Accept':'application/json','User-Agent':'COAD-source-review/1.0'})
        with urllib.request.urlopen(request,timeout=45) as response:data=response.read()
        obj=json.loads(data)
        path.write_bytes(data)
        registry[url]={'url':url,'retrieved_utc':datetime.now(timezone.utc).isoformat(),'sha256':hashlib.sha256(data).hexdigest(),'cache_file':path.name}
        reg_path.write_text(json.dumps(registry,indent=2)+'\n')
        time.sleep(.4)
        return obj
    path=ROOT/'results/COAD/03_PATIENT/20260919T114909Z_source_readiness_v2/gene_coverage.tsv'
    with path.open() as stream:rows=[r for r in csv.DictReader(stream,delimiter='\t') if r['exact_symbol_found']=='False']
    results=[]
    for row in rows:
        ident=row['human_gene_id'].split(':')[1]
        url='https://rest.genenames.org/fetch/entrez_id/'+ident
        data=fetch(url)
        docs=data['response']['docs']
        assert all(str(x.get('entrez_id'))==ident for x in docs)
        assert len(docs)<=1
        gene=docs[0] if docs else {}
        result={'human_gene_id':row['human_gene_id'],'original_symbol':row['gene'],
          'approved_symbol':gene.get('symbol'),'description':gene.get('name'),
          'hgnc_records':[
              {k:x.get(k) for k in ['hgnc_id','symbol','entrez_id','status','prev_symbol','alias_symbol','name','locus_type']} for x in docs],
          'source_urls':[url]}
        if ident=='102724197':
            result.update(ncbi_manual_review={'symbol':'LOC102724197','description':'inactive glutathione hydrolase 2',
              'method':'NCBI Gene page read via web retrieval; bulk ESummary access timed out',
              'source_url':'https://www.ncbi.nlm.nih.gov/gene/102724197','checked_date':'2026-09-19'})
            result['source_urls'].append('https://www.ncbi.nlm.nih.gov/gene/102724197')
        results.append(result)
        print(row['human_gene_id'],row['gene'],'->',gene.get('symbol'),'HGNC:',len(docs),flush=True)
    assert len(results)==22
    (CACHE/'gene_identity_records.json').write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    # These are the only previous symbols found in the author RNA row index.
    # Reverse lookup rejects a previous/alias symbol shared by multiple HGNC records.
    alias_checks=[]
    for original,old in [('GBA1','GBA'),('SLC35D4','TMEM241'),('SLC60A2','MFSD4B')]:
        target=next(r for r in results if r['original_symbol']==original)
        h=target['hgnc_records'][0]
        assert old in h['prev_symbol'] and h['status']=='Approved'
        query='(symbol:'+old+' OR prev_symbol:'+old+' OR alias_symbol:'+old+') AND status:Approved'
        url='https://rest.genenames.org/search/'+urllib.parse.quote(query)
        response=fetch(url)['response']
        assert response['numFound']==1 and response['docs'][0]['hgnc_id']==h['hgnc_id']
        alias_checks.append({'human_gene_id':target['human_gene_id'],'original_symbol':original,
          'rna_symbol':old,'hgnc_id':h['hgnc_id'],'status':'OFFICIAL_PREVIOUS_SYMBOL_UNIQUE',
          'source_url':url,'gene_source_url':target['source_urls'][0]})
    (CACHE/'alias_checks.json').write_text(json.dumps(alias_checks,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
