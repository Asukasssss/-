"""Local layout-only repair of untranslated source-state labels on report page5."""
from pathlib import Path
import json,shutil
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import pymupdf as fitz
R=Path(__file__).resolve().parents[2];O=R/'outputs/pdac-unified-report-v5-20260922-final'
c=json.loads((R/'configs/PDAC_report_v5.yaml').read_text(encoding='utf8'))
font='C:/Windows/Fonts/msyh.ttc';font_manager.fontManager.addfont(font)
plt.rcParams.update({'font.family':font_manager.FontProperties(fname=font).get_name(),'font.size':10,'pdf.fonttype':42,'svg.fonttype':'path','axes.unicode_minus':False,'axes.spines.top':False,'axes.spines.right':False})
df=pd.read_csv(O/'figure_data/05_mapping.tsv',sep='\t');left=df[df.unit.eq('入口代谢物')];right=df[~df.unit.eq('入口代谢物')]
fig,(ax,bx)=plt.subplots(1,2,figsize=(11.6,4.8));ax.barh([c['mapping_labels'].get(x,x) for x in left.label],left['count'],color='#087F8C')
for i,n in enumerate(left['count']):ax.text(n+.5,i,str(n),va='center')
ax.set_xlabel('入口代谢物项数');ax.set_xlim(0,34)
bx.bar(['直接关系','当前直接基因','条件或历史补充'],right['count'],color=['#087F8C','#337EAA','#9AA6B2']);bx.set_ylabel('数量（不同计数单位）');fig.tight_layout()
for ext in ['png','pdf','svg']:fig.savefig(O/'figures'/('05_mapping.'+ext),dpi=180,bbox_inches='tight')
plt.close(fig)
p=fitz.open(O/'PDAC_主报告.pdf');page=p[4];rect=fitz.Rect(30,94,812,438)
page.add_redact_annot(rect,fill=(1,1,1));page.apply_redactions();f=fitz.open(O/'figures/05_mapping.pdf');page.show_pdf_page(rect,f,0,keep_proportion=True);f.close()
box=fitz.Rect(30,451,410,546);page.add_redact_annot(box,fill=(1,1,1));page.apply_redactions()
assert page.insert_textbox(box,'结果与解释｜'+c['narratives']['mapping'][1],fontname='china-s',fontsize=11,color=(.12,.19,.25))>=0
for path in [O/'PDAC_主报告.md',O/'figures/05_mapping_caption_CN.md']:
 text=path.read_text(encoding='utf8').replace('7项无可靠关系','7项没有计划关系');path.write_text(text,encoding='utf8')
tmp=O/'PDAC_主报告_layout.tmp.pdf';p.save(tmp,garbage=4,deflate=True);p.close();tmp.replace(O/'PDAC_主报告.pdf')
p=fitz.open(O/'PDAC_主报告.pdf');p[4].get_pixmap(matrix=fitz.Matrix(1.3,1.3)).save(O/'rendered_pages/page_05.png');p.close()
shutil.copyfile(R/'configs/PDAC_report_v5.yaml',O/'repro/PDAC_report_v5.yaml')
(O/'analysis_spec.json').write_text(json.dumps(c,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
print('PAGE5_LABELS_TRANSLATED; statistics unchanged')
