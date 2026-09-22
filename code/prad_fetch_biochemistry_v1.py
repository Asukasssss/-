"""Fetch only public reviewed human UniProt biochemical entries for bounded mapping extension."""
import json, hashlib, time
from pathlib import Path
import requests

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/PRAD/02_MAPPING/20260922T142000Z_discovery_v1/sources/uniprot'
GENES='GSS GCLC GSR GPX1 GPX4 PGK1 PGK2 PGAM1 PGAM2 PHGDH AADAT AASS DHTKD1 GCK HK1 HK2 HK3 SLC2A1 SLC2A3 AKR7A2 AKR7A3 ADHFE1 PNP NT5C2 NT5C1A CTBS GNPNAT1 GLYCTK SHMT1 SHMT2 SDS SDSL GCAT GGT1 GGT5 GDPD1 GDPD3 ETNPPL SUGCT APRT BCAT1 BCAT2 NAGS ACY1 IMPA1 IMPA2 ISYNA1 NUDT5 NUDT9 ADPRH ADPRHL2 GALK2 STS SULT2A1 HSD3B7 CYP7B1 MGLL ABHD6 ABHD12 ENO1 ENO2 ENO3 PLA2G4A PLA2G6 CYP27A1 DPEP1 CNDP2'.split()

GENES += 'XDH SLC2A9 SLC22A12 MGAM GAA AMY1A AMY2A SULT1A1 HSD3B7'.split()

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    session=requests.Session();session.proxies={'http':'http://127.0.0.1:7890','https':'http://127.0.0.1:7890'}
    records=[]
    pending=[g for g in dict.fromkeys(GENES) if not (OUT/(g+'.json')).exists()]
    if (OUT/'fetch_manifest.json').exists():records=json.loads((OUT/'fetch_manifest.json').read_text())
    for start in range(0,len(pending),12):
        genes=pending[start:start+12]
        query='organism_id:9606 AND reviewed:true AND ('+' OR '.join('gene_exact:'+g for g in genes)+')'
        response=session.get('https://rest.uniprot.org/uniprotkb/search',params={'query':query,'format':'json','size':500},timeout=90);response.raise_for_status()
        entries=response.json()['results']
        assert 'rel="next"' not in response.headers.get('Link','')
        for g in genes:
            hits=[e for e in entries if any(v.get('geneName',{}).get('value')==g for v in e.get('genes',[]))]
            if len(hits)!=1:
                records.append(dict(gene=g,status='NEEDS_REVIEW',exact_primary_hits=len(hits)));continue
            p=OUT/(g+'.json');p.write_text(json.dumps(hits[0],ensure_ascii=False,indent=2),encoding='utf-8')
            records.append(dict(gene=g,status='DONE',accession=hits[0]['primaryAccession'],sha256=hashlib.sha256(p.read_bytes()).hexdigest(),release=response.headers.get('x-uniprot-release'),url=response.url))
            for c in hits[0].get('comments',[]):
                if c['commentType']=='CATALYTIC ACTIVITY':print(g,c['reaction']['name'])
        print('batch_complete',start+len(genes),flush=True)
    (OUT/'fetch_manifest.json').write_text(json.dumps(records,indent=2),encoding='utf-8')

if __name__=='__main__':main()
