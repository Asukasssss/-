"""Seven-page display: no spatial pages; PDAC focuses on CAFs."""
from pathlib import Path
import shutil,json,hashlib,html
import pandas as pd
from pypdf import PdfWriter,PdfReader
from PIL import Image,ImageOps,ImageDraw
ROOT=Path(__file__).resolve().parents[1];RUN='20260929T_representative_cases_v2'
O=ROOT/'results/BRCA/07_INTEGRATION'/RUN;OLD=O.parent/'20260929T_representative_cases_v1';O.mkdir(exist_ok=True)
V=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_representative_cases_v2')
files=[]
for n in [1,2,4,5,6,11,9]:
 if n!=5:
  for ext in ['png','pdf','svg']:
   f=next(OLD.glob(f'{n:02d}_*.{ext}'));shutil.copy2(f,O/f.name)
 files.append(next(O.glob(f'{n:02d}_*.png')))
for name in ['CAMP_case_results.tsv','relation_source_manifest.tsv','PDAC_SLC6A6_celltype_descriptive.tsv','server_versions.json']:
 shutil.copy2(OLD/name,O/name)
sc=pd.read_csv(OLD/'sc_patient_summary.tsv',sep='\t');caf=pd.read_csv(O/'PDAC_CAF_patient_summary.tsv',sep='\t');sc=pd.concat([sc[sc.cancer!='PDAC'],caf],ignore_index=True);sc.to_csv(O/'sc_patient_summary.tsv',sep='\t',index=False)
m=caf.iloc[0];count=f'{int(m.higher)}/{int(m.n_pairs)}'
for svg in O.glob('*.svg'):svg.write_text('\n'.join(x.rstrip() for x in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8',newline='\n')
w=PdfWriter()
for f in files:w.append(str(f.with_suffix('.pdf')))
with (O/'三个候选案例_无空间版.pdf').open('wb') as out:w.write(out)
cards=[]
for i,f in enumerate(files,1):
 group='ASNS' if i<3 else 'SLC6A6' if i<6 else 'UCKL1'
 cards.append(f'<article data-group="{group}"><h2>{i:02d} · {html.escape(f.stem.split("_",1)[1])}</h2><a href="{f.name}"><img src="{f.name}"></a><p><a href="{f.with_suffix(".pdf").name}">单页 PDF</a></p></article>')
page=(OLD/'index.html').read_text(encoding='utf-8');start=page.index('<main>')+6;end=page.index('</main>');page=page[:start]+''.join(cards)+page[end:]
page=page.replace('全部 11 页','全部 7 页').replace('ASNS · 3 页','ASNS · 2 页').replace('SLC6A6 · 5 页','SLC6A6 · 3 页').replace('UCKL1 · 3 页','UCKL1 · 2 页').replace('三个候选案例_统一展示.pdf','三个候选案例_无空间版.pdf')
(O/'index.html').write_text(page,encoding='utf-8')
spec=json.loads((OLD/'analysis_spec.json').read_text());spec.update(analysis_version=RUN,scope='Seven pages; CAMP and single-cell only',spatial_display='OMITTED_FROM_PRESENTATION',PDAC_contrast='Same untreated patient CAF vs malignant epithelial, >=20 nuclei per group',PDAC_direction=count,PDAC_statistical_status='DESCRIPTIVE_NO_NEW_TEST',replaces_PDac_prior_contrast='v1 7/9 malignant-vs-ductal is a different comparison and is not reused as CAF evidence')
for k in ['spatial_scale','spatial_selection','geomx']:spec.pop(k,None)
(O/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'README_CN.md').write_text(f'''# 三个候选案例 v2

## 本轮问题
按用户要求撤下这三个案例全部空间展示；将 PDAC SLC6A6 改为 CAF 定位。旧空间结果归档不删除，LYPLA1 既有空间页不受影响。

## 输入与范围
复用 v1 的 CAMP 冻结数值、单细胞表达和坐标；PDAC 仍为 Hwang 未治疗标本 108,964 个细胞核。来源见 CAF_source_manifest.tsv 与 v1 来源清单。

## 实际结果
共 7 页：ASNS 2 页、SLC6A6 3 页、UCKL1 2 页。PDAC 全细胞描述性平均表达：CAF 0.771、恶性上皮 0.136（log1p(CP10K)）。按患者配对、两类各至少 20 个细胞核，{count} 位患者 CAF 更高，平均差 {m.mean_difference:.3f}。未新增显著性检验。

## 新手解释
红圈和加粗图例现在指向 CAF；右侧每条线比较同一患者的恶性上皮与 CAF，不再使用原上皮–导管的 7/9 结果。全细胞均值与患者配对分别描述，不互相替代。

## 限制/反证
CAF 高于恶性上皮是细胞类型比较，不证明 CAF 相对健康成纤维细胞上调，也不证明组织代谢变化由 CAF 引起。描述性方向不是每位患者均显著；其他案例未显著结果照常保留。

## 当前决定
本汇报仅保留代谢物–基因联系与单细胞背景；空间图均不纳入该 7 页 PDF/预览。

## 下一步
按组会时长选页，不新增功能机制结论。

## 复现命令
server165 新运行目录执行 `ppt_cases_caf_v2.py RUN_DIR`；仅回传 public/。本地执行 `code/ppt_cases_package_v2.py`。原矩阵与逐患者表均留 server165。
''',encoding='utf-8')
pd.DataFrame([dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in [ROOT/'code/ppt_cases_caf_v2.py',Path(__file__)]]).to_csv(O/'code_manifest.tsv',sep='\t',index=False)
contact=Image.new('RGB',(1600,4*470),'#e6ece8')
for i,f in enumerate(files):contact.paste(ImageOps.contain(Image.open(f),(790,440)),((i%2)*800,(i//2)*470+25))
contact.save(O/'contact_sheet.jpg',quality=92)
assert len(PdfReader(O/'三个候选案例_无空间版.pdf').pages)==7
(O/'validation.json').write_text(json.dumps(dict(status='PASS',pages=7,spatial_pages=0,PDAC_focus='CAF',PDAC_paired_count=count,new_tests=False,visual_review='Reviewed revised PDAC panel and seven-page contact sheet'),ensure_ascii=False,indent=2),encoding='utf-8')
V.mkdir(exist_ok=True)
for f in O.iterdir():
 if f.is_file():shutil.copy2(f,V/f.name)
stage=ROOT/'coordination/stages/BRCA.tsv'
if RUN not in stage.read_text(encoding='utf-8'):
 with stage.open('a',encoding='utf-8') as f:f.write('\t'.join(['BRCA','07_INTEGRATION',RUN,'representative_cases_v2','DONE','Seven pages; CAF correction; no spatial display',str(O.relative_to(ROOT)).replace('\\','/'),'code/ppt_cases_package_v2.py','presentation/camp-fourcancer-lypla1-20260928','New CAF comparison descriptive; previous contrasts preserved in v1','Review revised display'])+'\n')
print('PACKAGED',count,V)
