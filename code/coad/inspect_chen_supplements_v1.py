"""Server-only DOCX structure audit. Export labels and dimensions, not patient rows."""
import hashlib,json,platform,zipfile,xml.etree.ElementTree as E
from pathlib import Path
import audit_relation39_sources_v1 as base
RUN='20260920T111614Z_supplement_audit_v1'
OUT=base.ROOT/'results/collaborative/COAD/B'/RUN
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
def txt(el):return ''.join(t.text or ''for t in el.findall('.//w:t',NS)).strip()
def main():
    assert Path.cwd()==OUT and (OUT/'.running').is_dir()
    inventory=json.loads((OUT/'supplement_inventory.json').read_text())
    assert len(inventory)==2 and all(x['status']=='DONE'for x in inventory)
    records=[]
    for x in inventory:
        p=base.DATA/x['filename']
        assert hashlib.sha256(p.read_bytes()).hexdigest()==x['sha256']
        with zipfile.ZipFile(p)as z:
            root=E.fromstring(z.read('word/document.xml'))
            # Only title and table header metadata; no body rows exported.
            titles=[]
            for para in root.findall('./w:body/w:p',NS):
                text=txt(para)
                if text.lower().startswith(('table s','supplementary table')):titles.append(text)
            headers=[]
            for table in root.findall('.//w:tbl',NS):
                first=table.find('./w:tr',NS)
                headers.append([txt(c)for c in first.findall('./w:tc',NS)])
            # Figure labels only; no figure text or numeric chart data exported.
            figure_labels=sorted(set(__import__('re').findall(r'(?:Fig(?:ure)?\.?\s*S\d+)',txt(root),__import__('re').I)))
            external=[]
            for name in z.namelist():
                if name.endswith('.rels'):
                    for item in E.fromstring(z.read(name)):
                        if item.get('TargetMode')=='External':external.append(dict(type=item.get('Type','').split('/')[-1],target=item.get('Target','')))
            records.append(dict(source_id=x['source_id'],table_titles=titles,table_headers=headers,figure_labels=figure_labels,external_links=external,embedded_objects=x['embedded_spreadsheet_files'],numeric_patient_rows_exported=False))
    (OUT/'structure_review.json').write_text(json.dumps(records,ensure_ascii=False,indent=2))
    meta=dict(run_id=RUN,status='DONE',scope='Two DOCX package structures and table headers inspected;not a numerical validation',python=platform.python_version(),source_code_commit='7f7363fb94ea95676e073b36ca8acc705d63a296',code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in [Path(__file__),OUT/'audit_relation39_sources_v1.py']},new_tests=False,patient_values_exported=False)
    (OUT/'execution_manifest.json').write_text(json.dumps(meta,indent=2))
    print(json.dumps(records,ensure_ascii=False))
if __name__=='__main__':main()
