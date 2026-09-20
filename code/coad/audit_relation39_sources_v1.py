"""Server-only supplementary audit; export inventory, never individual table rows."""
import hashlib,json,urllib.request,zipfile,xml.etree.ElementTree as E
from pathlib import Path
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
OUT=ROOT/'results/collaborative/COAD/B/20260920T105400Z_relation39_feasibility_v1'
DATA=ROOT/'data/candidates/coad_relation39_20260920'

def main():
    assert (OUT/'.running').exists();DATA.mkdir(exist_ok=True)
    rows=[];ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    for i in [1,2]:
        name=f'12943_2025_2359_MOESM{i}_ESM.docx'
        u='https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs12943-025-02359-x/MediaObjects/'+name
        p=DATA/name
        try:
            if not p.exists():
                with urllib.request.urlopen(u,timeout=60) as h:p.write_bytes(h.read())
            with zipfile.ZipFile(p)as z:
                assert z.testzip() is None
                x=E.fromstring(z.read('word/document.xml'))
                paragraphs=[''.join(t.itertext())for t in x.findall('.//w:p',ns)]
                # Only heading metadata and table dimensions. No individual rows or clinical values.
                labels=[]
                for para in x.findall('./w:body/w:p',ns):
                    t=''.join(v.text or '' for v in para.findall('.//w:t',ns)).strip()
                    if t.lower().startswith(('supplementary table','supplementary fig','table s','figure s')):labels.append(t[:180])
                tables=x.findall('.//w:tbl',ns)
                rows.append(dict(source_id='Chen2025_MOESM'+str(i),url=u,filename=name,status='DONE',bytes=p.stat().st_size,
                    sha256=hashlib.sha256(p.read_bytes()).hexdigest(),table_count=len(tables),table_rows=[len(t.findall('./w:tr',ns))for t in tables],
                    table_columns=[max([len(r.findall('./w:tc',ns))for r in t.findall('./w:tr',ns)]or[0])for t in tables],
                    headings=labels,embedded_spreadsheet_files=[n for n in z.namelist()if n.startswith('word/embeddings/')],
                    media_count=sum(n.startswith('word/media/')for n in z.namelist()),private_content_exported=False))
        except Exception as e:rows.append(dict(source_id='Chen2025_MOESM'+str(i),url=u,status='ACCESS_BLOCKED',reason=type(e).__name__+':'+str(e)))
    (OUT/'supplement_inventory.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows))
if __name__=='__main__':main()
