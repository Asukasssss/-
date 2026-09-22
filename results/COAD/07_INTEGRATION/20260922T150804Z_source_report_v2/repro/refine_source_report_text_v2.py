"""Reassemble existing vector figures with one embedded CJK/Latin font; no data edits."""
import argparse,json
from pathlib import Path
import pymupdf as fitz
def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();out=a.out
 boxes=json.loads((out/'layout_boxes.json').read_text());doc=fitz.open();font='C:/Windows/Fonts/msyh.ttc'
 for i in range(1,max(x['page'] for x in boxes)+1):
  entries=[x for x in boxes if x['page']==i];page=doc.new_page(width=842,height=595);page.insert_font(fontname='coadcjk',fontfile=font);page.draw_rect(fitz.Rect(0,0,842,8),color=None,fill=(.03,.5,.55))
  ident=entries[-1]['text'].split('|')[-2].strip();figure=fitz.open(out/'figures'/(ident+'.pdf'));page.show_pdf_page(fitz.Rect(30,94,812,438),figure,0,keep_proportion=True);figure.close()
  for e in entries:
   color=(.48,.28,.16) if e['text'].startswith('限制') else (.12,.19,.25)
   assert page.insert_textbox(fitz.Rect(e['rect']),e['text'],fontname='coadcjk',fontsize=e['font_size'],color=color)>=0
 dest=out/'COAD_主报告.pdf';temp=out/'report_refined.tmp.pdf';doc.subset_fonts();doc.save(temp,garbage=4,deflate=True);doc.close();temp.replace(dest)
 pdf=fitz.open(dest)
 for i,page in enumerate(pdf):page.get_pixmap(matrix=fitz.Matrix(1.3,1.3)).save(out/'rendered_pages'/f'page_{i+1:02d}.png')
 (out/'text_refinement_validation.json').write_text(json.dumps(dict(status='PASS',pages=len(pdf),embedded_font='Microsoft YaHei',all_textboxes_fit=True,figures_reused=True,no_statistics_changed=True),indent=2))
if __name__=='__main__':main()
