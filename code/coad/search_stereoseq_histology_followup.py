import urllib.request,json,re,concurrent.futures
from pathlib import Path
def get(u):
    with urllib.request.urlopen(u,timeout=30) as r:return r.read().decode('utf-8')
def detail(i):
    u=f'https://db.cngb.org/stomics/tissue_section/STTS{i:07d}'
    try:
        t=get(u);d=json.JSONDecoder().raw_decode(t.split('window.__INITIAL_STATE__=')[1])[0]['PublicPages']['tissueSectionView']
        return {'url':u,'data':d,'image_tags':[x[:200] for x in re.findall('<img[^>]+>',t) if 'logo' not in x and 'qr-code' not in x]}
    except Exception as e:return {'url':u,'error':str(e)}
with concurrent.futures.ThreadPoolExecutor(4) as ex: rows=list(ex.map(detail,range(650,666)))
Path('runtime/histology_details_all.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
for r in rows:print(r['url'],r.get('data',{}).get('title'),r.get('data',{}).get('tissue_section_processing',{}).get('staining_protocol'),r.get('image_tags'),r.get('error',''))
for u in ['https://ftp.cngb.org/pub/CNSA/data4/CNP0002432/','https://api.github.com/repos/STOmics/StereoMMv1/releases','https://api.github.com/repos/STOmics/StereoMMv1/branches','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12100622/fullTextXML']:
    try:
        t=get(u)
        if 'ftp.cngb' in u:print('FTP_LINKS',re.findall(r'href="([^"]+)"',t))
        elif 'github' in u:print(u,t[:4500])
        else:
            Path('runtime/stereomm.xml').write_text(t,encoding='utf-8');print('PAPER_URLS', sorted(set(re.findall(r'https?[^<>"\s]+',t))))
    except Exception as e:print(u,type(e).__name__,str(e))
