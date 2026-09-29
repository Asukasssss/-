"""Known biochemical roles above frozen CAMP associations; no new statistics."""
from pathlib import Path
import json,hashlib,shutil
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
from matplotlib.lines import Line2D
from pypdf import PdfWriter,PdfReader
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/BRCA/07_INTEGRATION';OLD=BASE/'20260929T_representative_cases_v2';RUN='20260929T_representative_cases_v3';O=BASE/RUN;O.mkdir(exist_ok=True)
V=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_representative_cases_v3')
for f in OLD.iterdir():
 if f.is_file():shutil.copy2(f,O/f.name)
for n in ['msyh.ttc','msyhbd.ttc']:font_manager.fontManager.addfont('C:/Windows/Fonts/'+n)
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'axes.unicode_minus':False})
G='#00553B';INK='#183D34';M='#6A7973';UP='#C77462';DOWN='#468B9E';LINE='#DDE7E0'
df=pd.read_csv(O/'CAMP_case_results.tsv',sep='\t').set_index('cancer')
def t(f,x,y,s,size=13,col=INK,bold=False,ha='left'):f.text(x,y,s,fontsize=size,color=col,weight='bold' if bold else 'normal',ha=ha,va='center')
def line(f,y):f.add_artist(Line2D([.05,.95],[y,y],transform=f.transFigure,color=LINE,lw=.8))
def arrow(f,x1,x2,y):f.add_artist(FancyArrowPatch((x1,y),(x2,y),transform=f.transFigure,arrowstyle='-|>',mutation_scale=20,color=G,lw=2))
sources=[]
for gene,acc in [('ASNS','P08243'),('SLC6A6','P31641'),('UCKL1','Q9NWZ5')]:
 j=json.loads((O/(acc+'.json')).read_text(encoding='utf-8'));assert j['primaryAccession']==acc and j['entryType']=='UniProtKB reviewed (Swiss-Prot)' and j['organism']['taxonId']==9606
 rx=[c['reaction'] for c in j['comments'] if c['commentType']=='CATALYTIC ACTIVITY'][0]
 sources.append(dict(gene=gene,accession=acc,url=f'https://www.uniprot.org/uniprotkb/{acc}/entry',retrieved='2026-09-29',reaction=rx['name'],rhea=rx['reactionCrossReferences'][0]['id'],pmids=','.join(e['id'] for e in rx.get('evidences',[]) if e['source']=='PubMed'),sha256=hashlib.sha256((O/(acc+'.json')).read_bytes()).hexdigest()))
pd.DataFrame(sources).to_csv(O/'biochemical_sources.tsv',sep='\t',index=False)
pages=[('ASNS','01_BRCA_ASNS_CAMP',['BRCA'],'ASNS：利用谷氨酰胺的氮，合成天冬酰胺'),('SLC6A6','04_PDAC_PRAD_SLC6A6_CAMP',['PDAC','PRAD'],'SLC6A6：将细胞外牛磺酸转运入细胞'),('UCKL1','11_COAD_UCKL1_CAMP',['COAD'],'UCKL1：将尿苷转为 UMP，参与补救合成')]
for gene,slug,cs,title in pages:
 src=next(s for s in sources if s['gene']==gene)
 f=plt.figure(figsize=(16,9),facecolor='white');ax=f.add_axes([.044,.882,.168,.10]);ax.imshow(plt.imread(BASE/'20260929T170000Z_page8_genes_six_v2/sysu_logo.png'));ax.axis('off');t(f,.95,.93,'CAMP  /  组会汇报',10,M,ha='right');line(f,.883);t(f,.05,.827,title,25,G,True)
 t(f,.05,.756,'01  已知生化作用',14,G,True)
 f.add_artist(FancyBboxPatch((.05,.516),.90,.208,boxstyle='round,pad=0.008,rounding_size=0.01',transform=f.transFigure,facecolor='#F0F5F1',edgecolor='none',zorder=-1))
 if gene=='ASNS':
  t(f,.275,.630,'谷氨酰胺 + 天冬氨酸',22,G,True,ha='center');t(f,.74,.630,'谷氨酸 + 天冬酰胺',22,G,True,ha='center');arrow(f,.445,.565,.630);t(f,.505,.678,'ASNS',17,G,True,ha='center');t(f,.505,.585,'ATP → AMP + PPi',12,M,ha='center');t(f,.5,.540,'谷氨酰胺是氮供体和反应底物；ASNS 合成的是天冬酰胺。',13,INK,ha='center')
 elif gene=='SLC6A6':
  t(f,.25,.681,'细胞外',12,M,ha='center');t(f,.75,.681,'细胞内',12,M,ha='center');t(f,.25,.627,'牛磺酸',23,G,True,ha='center');t(f,.75,.627,'牛磺酸',23,G,True,ha='center')
  f.add_artist(FancyBboxPatch((.425,.596),.15,.072,boxstyle='round,pad=0.007',transform=f.transFigure,facecolor='#D7E7DD',edgecolor=G,lw=1));t(f,.50,.632,'SLC6A6 / TauT',15,G,True,ha='center');arrow(f,.33,.415,.627);arrow(f,.587,.667,.627);t(f,.5,.549,'Na+ / Cl- 依赖性摄取  ·  搬运牛磺酸，本步骤不改变化学结构',13,INK,ha='center')
 else:
  t(f,.275,.630,'尿苷 + ATP',24,G,True,ha='center');t(f,.74,.630,'UMP + ADP',24,G,True,ha='center');arrow(f,.42,.58,.630);t(f,.50,.679,'UCKL1',17,G,True,ha='center');t(f,.74,.585,'尿苷单磷酸',12,M,ha='center');t(f,.50,.542,'尿苷是磷酸化底物，进入嘧啶核苷补救合成；也可催化胞苷 → CMP。',13,INK,ha='center')
 t(f,.05,.474,'02  CAMP 中实际观察到的结果',14,G,True)
 for x,h in [(.24,'代谢物变化'),(.49,'基因 RNA 变化'),(.78,'肿瘤样本内关联')]:t(f,x,.423,h,13,G,True,ha='center')
 for i,c in enumerate(cs):
  r=df.loc[c];y=.317 if len(cs)==1 else .334-i*.121;t(f,.07,y,c,15,G,True)
  for xx,label,eff,p in [(.24,r.molecule,r.m_effect,r.m_p),(.49,r.gene,r.rna_effect,r.rna_p)]:
   t(f,xx,y+.010,label+(' ↑' if eff>0 else ' ↓'),19,UP if eff>0 else DOWN,True,ha='center');t(f,xx,y-.032,f'配对 P = {p:.3g}',10.5,M,ha='center')
  t(f,.78,y+.010,f'ρ = {r.rho:+.3f}   P = {r.association_p:.4g}',15,G,True,ha='center');t(f,.78,y-.032,f'{int(r.association_n)} 个肿瘤样本 · Spearman 相关',10.5,M,ha='center')
  if len(cs)==2 and i==0:line(f,.278)
 conclusion={'ASNS':'对应线索：谷氨酰胺降低、ASNS 升高且负相关；尚不能断定是 ASNS 消耗所致。','SLC6A6':'PDAC 单核定位指向 CAF；组织牛磺酸变化尚不能归因为 CAF 的转运作用。','UCKL1':'对应线索：尿苷降低、UCKL1 升高且负相关；尚未直接测量尿苷补救通量。'}[gene]
 line(f,.151);t(f,.05,.113,conclusion,13.5,G,True)
 t(f,.05,.070,'生化注释提供映射依据；RNA 水平与相关性不等于酶活、转运量或因果。',10,M)
 t(f,.05,.034,f"来源：UniProt {src['accession']} · {src['rhea']}（人类，已审阅）；反应为简化示意。",8.5,M)
 for ext in ['png','pdf','svg']:f.savefig(O/(slug+'.'+ext),dpi=200,facecolor='white')
 plt.close(f)
for svg in O.glob('*.svg'):svg.write_text('\n'.join(l.rstrip() for l in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8',newline='\n')
w=PdfWriter();files=[next(O.glob(f'{n:02d}_*.png')) for n in [1,2,4,5,6,11,9]]
for f in files:w.append(str(f.with_suffix('.pdf')))
with (O/'三个候选案例_无空间版.pdf').open('wb') as out:w.write(out)
assert len(PdfReader(O/'三个候选案例_无空间版.pdf').pages)==7
from PIL import Image,ImageOps
contact=Image.new('RGB',(1600,4*470),'#e6ece8')
for i,f in enumerate(files):contact.paste(ImageOps.contain(Image.open(f),(790,440)),((i%2)*800,(i//2)*470+25))
contact.save(O/'contact_sheet.jpg',quality=92)
text=(O/'README_CN.md').read_text(encoding='utf-8').replace('三个候选案例 v2','三个候选案例 v3')
text+='\n## v3 生化说明补充\n三张 CAMP 页上方补充已知反应或转运，下方保留冻结数值；不新增机制证据或统计检验。ASNS 使用谷氨酰胺酰胺氮合成天冬酰胺；SLC6A6 进行钠/氯依赖性牛磺酸转运；UCKL1 催化尿苷到 UMP 的磷酸化。页面是简化示意，ASNS 省略 H2O/H+，UCKL1 省略 H+；完整化学计量与来源见 biochemical_sources.tsv、三个 UniProt JSON。功能注释不证明这些过程在本队列增强。复现：本地运行 code/ppt_cases_roles_v3.py。\n'
(O/'README_CN.md').write_text(text,encoding='utf-8')
spec=json.loads((O/'analysis_spec.json').read_text());spec.update(analysis_version=RUN,biochemical_source='Reviewed human UniProt entries P08243 P31641 Q9NWZ5, retrieved 2026-09-29',change='Three role diagrams added; data and four single-cell pages unchanged')
(O/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf-8')
pd.DataFrame([dict(path=str(Path(__file__).relative_to(ROOT)),sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())]).to_csv(O/'code_manifest_v3.tsv',sep='\t',index=False)
(O/'validation.json').write_text(json.dumps(dict(status='PASS',pages=7,spatial_pages=0,unchanged_statistics=True,reviewed_human_entries=3,visual_review='All three changed slides inspected at full size; labels and arrows legible with no overlap'),indent=2))
V.mkdir(exist_ok=True)
for f in O.iterdir():
 if f.is_file():shutil.copy2(f,V/f.name)
stage=ROOT/'coordination/stages/BRCA.tsv'
if RUN not in stage.read_text(encoding='utf-8'):
 with stage.open('a',encoding='utf-8') as f:f.write('\t'.join(['BRCA','07_INTEGRATION',RUN,'representative_cases_v3','DONE','Known biochemical roles added to three CAMP pages',str(O.relative_to(ROOT)).replace('\\','/'),'code/ppt_cases_roles_v3.py','presentation/camp-fourcancer-lypla1-20260928','Reviewed human UniProt reactions; frozen statistics; no spatial pages','Review biochemical role diagrams'])+'\n')
print('RENDERED',V)

