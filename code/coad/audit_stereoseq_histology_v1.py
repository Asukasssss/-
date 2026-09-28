"""Audit public supplementary assets on server165; no patient values are exported."""
import argparse, hashlib, json, zipfile, io
from pathlib import Path
from xml.etree import ElementTree as ET
from PIL import Image, ImageDraw

def main():
    p=argparse.ArgumentParser(); p.add_argument('--source',required=True); p.add_argument('--out',required=True)
    a=p.parse_args(); src=Path(a.source); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    records=[]
    for f in sorted(src.iterdir()):
        if f.is_file():
            records.append(dict(file=f.name,bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
    z=zipfile.ZipFile(src/'PMC11599708_supplementary.zip'); assert z.testzip() is None
    with zipfile.ZipFile(src/'41467_2024_54710_MOESM9_ESM.xlsx') as w:
        root=ET.fromstring(w.read('xl/workbook.xml'))
        sheets=[e.attrib['name'] for e in root.iter() if e.tag.endswith('}sheet')]
        media=[n for n in w.namelist() if n.startswith('xl/media/')]
        canvas=Image.new('RGB',(1200,240*((len(media)+3)//4)),'white'); draw=ImageDraw.Draw(canvas)
        for i,n in enumerate(media):
            im=Image.open(io.BytesIO(w.read(n))).convert('RGB'); im.thumbnail((290,210))
            x=(i%4)*300;y=(i//4)*240;canvas.paste(im,(x,y+25));draw.text((x+5,y+5),n,fill='black')
        canvas.save(out/'source_data_media_contact.png')
    report=dict(source_files=records,archive_integrity='PASS',source_data_sheets=sheets,source_data_media=media,
                note='Media contact sheet is for visual source classification; not a new tissue analysis.')
    (out/'asset_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
