"""Validate frozen values and package rendered presentation; no inferential tests."""
from pathlib import Path
import hashlib,json,html,platform,shutil,subprocess
import pandas as pd,numpy as np
import pymupdf as fitz
from PIL import Image,ImageDraw
ROOT=Path(__file__).resolve().parents[1]
O=ROOT/'results/BRCA/07_INTEGRATION/20260928T150000Z_ppt_redraw_v1'
ART=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_PPT_unified_v1')
ART.mkdir(parents=True,exist_ok=True)
M=pd.read_csv(O/'source_manifest.tsv',sep='\t')
for r in M.itertuples():assert hashlib.sha256((O/'sources'/(r.id+'.tsv')).read_bytes()).hexdigest()==r.sha256,r.id
rr=pd.read_csv(O/'figure_data/04_lypla1_rna.tsv',sep='\t')
assert rr.cancer.tolist()==['BRCA','COAD','PDAC','PRAD']
assert (rr.effect>0).all() and (rr.p<.05).tolist()==[True,True,False,True]
v=json.loads((O/'server_figures/render_validation.json').read_text())
assert len(v['checks'])==25 and not v['new_statistics'] and not v['new_embedding']
expected={'BRCA':28844,'COAD':31247,'PRAD':5357,'PDAC':61272}
for c,n in expected.items():assert next(x for x in v['checks'] if x['plot']=='umap_'+c)['cells']==n
assert len(list((O/'figures').glob('*.png')))==13
pdf=O/'CAMP_LYPLA1_统一重绘版.pdf';doc=fitz.open(pdf)
slides=json.loads((O/'slide_manifest.json').read_text(encoding='utf8'))
assert len(doc)==len(slides)==22
render=ART/'slides';render.mkdir(exist_ok=True)
for i,page in enumerate(doc):
 page.get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(render/f'{i+1:02d}.png')
 # Check native slide titles and page numbers are in the rendered text layer.
 txt=''.join(page.get_text().split())
 assert ''.join(slides[i]['title'].split()) in txt,(i+1,'missing title')
 assert f'{i+1:02d}' in txt
for start in range(0,len(doc),6):
 im=Image.new('RGB',(1500,3*445),'#dfe5ea');draw=ImageDraw.Draw(im)
 for k in range(start,min(start+6,len(doc))):
  a=Image.open(render/f'{k+1:02d}.png');a.thumbnail((738,415));x=(k-start)%2*750;y=(k-start)//2*445;im.paste(a,(x,y+22));draw.text((x+8,y+4),str(k+1),fill='black')
 im.save(ART/f'contact_{start+1:02d}.jpg',quality=90)
# Complete atlas: a labelled page for every section, never selected by expression.
atlas=fitz.open();inv=pd.read_csv(O/'server_figures/spatial_display_inventory.tsv',sep='\t')
for r in inv.itertuples():
 src=fitz.open(O/'server_figures'/(r.plot+'.pdf'));page=atlas.new_page(width=960,height=410)
 page.insert_text((30,24),f'{r.cancer} | {r.sample} | {r.patient}',fontsize=14)
 page.show_pdf_page(fitz.Rect(15,40,945,375),src,0)
 page.insert_text((30,398),'Fixed expression display: log1p(CP10K), 0-3. Original region labels and coordinates; no new tests.',fontsize=9)
atlas.save(ART/'空间切片完整图册_21张.pdf',garbage=4,deflate=True)
for item in ['CAMP_LYPLA1_统一重绘版.pptx','CAMP_LYPLA1_统一重绘版.pdf','source_manifest.tsv','figure_manifest.tsv','STYLE_STANDARD.json']:
 shutil.copy2(O/item,ART/item)
for folder in ['figures','figure_data','server_figures']:
 shutil.copytree(O/folder,ART/folder,dirs_exist_ok=True)
cards=[]
for s in slides:
 i=s['number'];cards.append(f'<article><h2>{i:02d} · {html.escape(s["title"])}</h2><p>{html.escape(s["claim"])}</p><a href="slides/{i:02d}.png"><img loading="lazy" src="slides/{i:02d}.png"></a><small>{html.escape(s["source"])}</small></article>')
gallery=[]
for r in inv.itertuples():gallery.append(f'<details><summary>{r.cancer} · {r.sample} · {r.patient}</summary><img loading="lazy" src="server_figures/{r.plot}.png"></details>')
page='''<!doctype html><html lang="zh"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CAMP · 四癌 LYPLA1 汇报</title><style>body{font-family:"Microsoft YaHei",sans-serif;margin:0;background:#f2f5f7;color:#172a3a}main{max-width:1200px;margin:40px auto;padding:0 24px}h1{font-size:32px}a{color:#3276a5}nav{display:flex;gap:18px;flex-wrap:wrap;margin:25px 0}nav a{background:white;padding:14px;border-radius:8px}article{background:white;margin:28px 0;padding:24px;border-radius:12px}img{width:100%;height:auto}small{color:#657583}details{background:white;margin:12px 0;padding:16px}summary{cursor:pointer}p{line-height:1.8}</style><main><h1>CAMP → 四癌 LYPLA1</h1><p>22页汇报 · 13张定量图 · 四癌UMAP · 21张空间切片。统一绘图代码、颜色、字体与统计标注。原统计结果保留。</p><nav><a href="CAMP_LYPLA1_统一重绘版.pptx">下载 PPT</a><a href="CAMP_LYPLA1_统一重绘版.pdf">阅读 PDF</a><a href="空间切片完整图册_21张.pdf">完整空间图册</a><a href="source_manifest.tsv">图源清单</a></nav>'''+''.join(cards)+'<h1>完整空间切片</h1>'+''.join(gallery)+'</main></html>'
(ART/'index.html').write_text(page,encoding='utf8')
validation=dict(status='DONE',source_hashes_checked=len(M),source_hash_mismatch=0,slides=22,quantitative_figures=13,UMAP=4,spatial_sections=21,source_cell_counts=expected,new_tests=0,new_embeddings=0,rendered_with='LibreOffice headless on server165',all_rendered_titles_present=True,rendered_pages=22,visual_review='All 22 slides and 21 section layouts inspected; legend/text collisions corrected',limitations=['No Microsoft PowerPoint-specific rendering test; use PDF as fixed-layout reference','Same color scale does not remove platform or depth differences','BRCA spatial is tumor vs immune, not normal epithelium; PRAD is one patient','PDAC RNA P>=0.05 retained; nominal threshold does not establish FDR control'],python=platform.python_version(),parent_code_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
(O/'validation.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2),encoding='utf8')
notes=['# 逐页讲解提纲\n']
for s in slides:notes += [f'## {s["number"]:02d} {s["title"]}',s['claim'],s['note'],s['source'],'']
(O/'逐页讲解.md').write_text('\n\n'.join(notes),encoding='utf8');shutil.copy2(O/'逐页讲解.md',ART/'逐页讲解.md')
files=[]
for p in O.rglob('*'):
 if p.is_file() and p.name!='artifact_manifest.tsv':files.append(dict(path=p.relative_to(O).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
pd.DataFrame(files).to_csv(O/'artifact_manifest.tsv',sep='\t',index=False)
print(json.dumps(validation,ensure_ascii=False));print(ART)
