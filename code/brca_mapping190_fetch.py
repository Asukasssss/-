"""Fetch public human protein annotation only; no patient data."""
import concurrent.futures, json, sys, time
from pathlib import Path
import requests

GENES='SDHA SDHB SDHC SDHD SUCLG1 SUCLG2 SUCLA2 SLC13A3 SLC25A10 MSRA ACY1 SLC22A2 SLC47A1 SLC47A2 CRAT SLC25A20 CPT1A CPT1B CPT2 ACSL1 ACSL3 ACSL4 ACSL5 ACSL6 ACSM1 ACSM2A ACSM2B ACSM3 LPCAT1 LPCAT2 LPCAT3 LPCAT4 LYPLA1 LYPLA2 ENPP2 MGLL ABHD6 ABHD12 GGT1 GGCT D2HGDH L2HGDH ALOX15 HPGD CYP27A1 HSD3B7 CYP7B1 SLC5A9 SLC5A10 NT5C2 PNP SORD SLC22A4 GOT1 GOT2 LDHA LDHB PGP ALDH1A1 AKR1D1'.split()

def fetch(gene, folder):
    target=folder/(gene+'.json')
    if target.exists(): return gene, 'CACHED'
    query=f'(gene_exact:{gene}) AND (organism_id:9606) AND (reviewed:true)'
    for attempt in range(3):
        try:
            r=requests.get('https://rest.uniprot.org/uniprotkb/search',params={'query':query,'format':'json','size':50},timeout=40)
            r.raise_for_status(); data=r.json()
            assert data.get('results'), gene
            target.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
            return gene, 'OK'
        except Exception as e:
            if attempt==2: return gene, str(e)
            time.sleep(2)

if __name__=='__main__':
    folder=Path(sys.argv[1]); folder.mkdir(parents=True,exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:
        for result in ex.map(lambda g:fetch(g,folder),GENES): print(*result,flush=True)
