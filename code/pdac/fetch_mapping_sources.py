"""Fetch public biochemical annotations; never reads patient data."""
import csv,hashlib,json,time,urllib.request,urllib.parse
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[2]
CACHE=ROOT/'runtime/pdac_mapping_v1';CACHE.mkdir(parents=True,exist_ok=True)
REG=CACHE/'sources.json'
sources=json.loads(REG.read_text()) if REG.exists() else {}
def fetch(url,name):
    p=CACHE/name
    if p.exists() and url in sources and hashlib.sha256(p.read_bytes()).hexdigest()==sources[url]['sha256']:
        return p.read_text(encoding='utf-8')
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'PDAC-biochemical-review/1.0'}),timeout=90) as f:data=f.read()
            p.write_bytes(data);sources[url]={'path':name,'url':url,'sha256':hashlib.sha256(data).hexdigest(),'retrieved_utc':datetime.now(timezone.utc).isoformat()}
            REG.write_text(json.dumps(sources,indent=2));print('FETCHED',name,len(data),flush=True);return data.decode()
        except Exception:
            if attempt==2:raise
            time.sleep(2)
def main():
    rows=list(csv.DictReader((ROOT/'results/PDAC/01_CAMP/20260920T051000Z_frozen_v1/significant_features.tsv').open(),delimiter='\t'))
    ids=sorted({r['kegg_id'] for r in rows if r['kegg_id']})
    for start in range(0,len(ids),8):
        fetch('https://rest.kegg.jp/get/'+'+'.join(ids[start:start+8]),f'kegg_{start}.txt')
    for name in ['chebi_pH7_3_mapping.tsv','chebiId_name.tsv']:
        fetch('https://ftp.expasy.org/databases/rhea/tsv/'+name,name)
    query={'query':'organism_id:9606 AND reviewed:true','format':'json','fields':'accession,gene_names,organism_id,xref_geneid,cc_catalytic_activity,cc_subunit,cc_subcellular_location,cc_function'}
    fetch('https://rest.uniprot.org/uniprotkb/stream?'+urllib.parse.urlencode(query),'human_reviewed.json')
if __name__=='__main__':main()
