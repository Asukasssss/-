"""Package public representative-case figures without copying private matrices."""
from pathlib import Path
import json,hashlib,html,shutil,subprocess
import pandas as pd
from PIL import Image,ImageOps,ImageDraw
from pypdf import PdfWriter,PdfReader
ROOT=Path(__file__).resolve().parents[1]
RUN='20260929T_representative_cases_v1'
O=ROOT/'results/BRCA/07_INTEGRATION'/RUN
V=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_representative_cases_v1')
order=[1,2,3,4,5,6,7,8,11,9,10]
files=[next(O.glob(f'{n:02d}_*.png')) for n in order]
assert len(files)==11
for svg in O.glob('*.svg'):
 svg.write_text('\n'.join(s.rstrip() for s in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8',newline='\n')
writer=PdfWriter()
for f in files:writer.append(str(f.with_suffix('.pdf')))
with (O/'三个候选案例_统一展示.pdf').open('wb') as out:writer.write(out)
assert len(PdfReader(O/'三个候选案例_统一展示.pdf').pages)==11
rows=[]
for p in sorted((ROOT/'code').glob('ppt_cases_*.py')):
 rows.append(dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
pd.DataFrame(rows).to_csv(O/'code_manifest.tsv',sep='\t',index=False)
spec=dict(analysis_version=RUN,status='DONE',scope='Three representative contexts; not a formal test of cancer specificity',excluded_presentation='FUSCC and Tang ASNS external association results',new_inferential_tests=False,p_values='Frozen CAMP, COAD sign test and PRAD paired test only; ASNS and PDAC single-cell descriptive',statistical_unit='Patient; spatial spots descriptive only',minimum_cells_per_group_per_patient=20,expression='log1p(gene raw count / full gene-expression library * 10000)',umap='Original accepted coordinates and author annotations reused; no new embedding or batch correction',umap_scale='0 to each panel positive-expression 99th percentile; zeros light gray',spatial_scale='0 to 3 log1p(CP10K), viridis, zero light gray; saturated values indicated',spatial_selection='BRCA and PRAD same user-selected LYPLA1 sections; COAD first two colon/cecum source cases without using target expression',geomx='Frozen LOQ rule (>10 ROIs above global negative-probe LOQ), Q3 primary; CPM and author VST descriptive sensitivity; six matched patients',seed=20260929,private_data='All per-cell, per-patient and per-ROI tables remain server165',source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
(O/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2),encoding='utf-8')
readme='''# 三个候选案例：统一展示

## 本轮问题
展示 BRCA 谷氨酰胺–ASNS、PDAC/PRAD 牛磺酸–SLC6A6、COAD 尿苷–UCKL1，保持 LYPLA1 的展示风格。按用户要求，不展示 FUSCC/Tang 外部关联。

## 输入与范围
CAMP 原有冻结结果；Wu、Uhlitz、Hwang、PRAD24 的真实全细胞/单核表达和既有坐标；CTA、Erickson、Valdeolivas 病理标注切片；Bell GeoMx 区域数据。所有原矩阵、逐患者/细胞/ROI 数值留在 server165。本包只含代码、汇总与图。

## 实际结果
|案例|CAMP 肿瘤标本相关|同患者细胞比较|空间补充|
|---|---|---|---|
|BRCA ASNS|谷氨酰胺下降、ASNS 上升；ρ=-0.440，P=0.0005|5/8 恶性上皮更高，描述性|两位患者 CTA 切片|
|PDAC SLC6A6|牛磺酸与 RNA 上升；ρ=0.494，P=0.0275|7/9 恶性上皮更高，描述性|GeoMx 4/6 患者癌上皮高于正常导管；三种尺度方向计数一致|
|PRAD SLC6A6|牛磺酸与 RNA 下降；ρ=0.141，P=0.1777|髓系肿瘤–癌旁 4/11 更高，P=0.6377|同一患者两张病理切片|
|COAD UCKL1|尿苷下降、UCKL1 上升；ρ=-0.491，P=0.0041|7/9 CNA 上皮更高，既有精确符号检验 P=0.1797|两病例；第二例保留非肿瘤上皮参照|

## 新手解释
每个案例先说明代谢物–基因联系，再用同一套全细胞 UMAP 展示表达和细胞身份，最后展示空间位置。每条患者连线是一位患者两组细胞的平均表达，不把细胞数量当作患者重复。相关系数不是倍数；RNA 不等于酶活或代谢通量。

## 限制/反证
三例是代表性背景，不是经癌种间交互检验证实的癌特异机制。组织相关不证明因果；单细胞的上皮参照与髓系参照分别标明，不可混称正常上皮。BRCA ASNS 与 PDAC SLC6A6 本轮只新增描述性患者汇总，未新增 P 值。COAD 的冻结 P 是符号检验，并非对图中均值做 t 检验。PRAD 单细胞阴性如实保留。

空间切片使用全库计数标准化、原坐标和作者病理标注，无平滑或插补，也未新增 spot 级显著性检验。BRCA 未标注区不等于正常上皮；PRAD 两张切片不是两位患者；COAD S1 缺合格非肿瘤上皮参照。GeoMx 是选定 ROI，不能伪装成完整切片热图。每张单细胞图色标上限不同且明确显示，不能靠颜色跨癌比较绝对丰度。

## 当前决定
11 页全部保留供组会选用：ASNS 3 页、SLC6A6 5 页、UCKL1 3 页。采用同一字体、颜色含义、右侧 UMAP 图例、患者连线及 H&E–病理区域–表达三联结构。

## 下一步
根据汇报时长选择主讲与备份页；当前包不包含新的功能机制验证。

## 复现命令
本地运行 `code/ppt_cases_relations_v1.py`。在 server165 新运行目录中依次运行 `ppt_cases_sc_prepare_v1.py RUN_DIR`、`ppt_cases_sc_render_v1.py RUN_DIR`（同目录需有仓库 `ppt_four_sc_render_v3.py`）、`ppt_cases_spatial_v1.py RUN_DIR`。仅回传 `public/`，再本地运行 `code/ppt_cases_package_v1.py`。脚本沿用明确路径、原阈值、标本标识和固定种子；逐细胞/ROI缓存不回传。
'''
(O/'README_CN.md').write_text(readme,encoding='utf-8')
cards=[]
for i,f in enumerate(files,1):
 name=f.stem.split('_',1)[1];group='ASNS' if i<=3 else 'SLC6A6' if i<=8 else 'UCKL1'
 cards.append(f'<article id="p{i}" data-group="{group}"><h2>{i:02d} · {html.escape(name)}</h2><a href="{f.name}" target="_blank"><img src="{f.name}" alt="{html.escape(name)}"></a><p><a href="{f.with_suffix(".pdf").name}">单页 PDF</a> · <a href="{f.with_suffix(".svg").name}">SVG</a></p></article>')
page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CAMP · 三个候选案例</title><style>body{margin:0;background:#eef3f0;color:#174338;font:16px "Microsoft YaHei",sans-serif}header{position:sticky;top:0;background:#fffffff5;padding:16px 4vw;z-index:2;border-bottom:1px solid #cddbd4;display:flex;align-items:center;gap:16px;flex-wrap:wrap}h1{font-size:20px;margin:0 auto 0 0}button,a{color:#00553b}button{background:white;border:1px solid #bbcfc4;border-radius:6px;padding:8px 14px;cursor:pointer}button.active{background:#00553b;color:white}main{max-width:1600px;margin:24px auto;padding:0 20px}article{margin-bottom:28px;background:white;border-radius:12px;overflow:hidden;border:1px solid #d8e3dc}h2{font-size:17px;padding:0 22px}img{display:block;width:100%;height:auto}p{padding:0 22px}article[hidden]{display:none}</style><header><h1>CAMP · 三个候选案例</h1><button class="active" data-filter="all">全部 11 页</button><button data-filter="ASNS">ASNS · 3 页</button><button data-filter="SLC6A6">SLC6A6 · 5 页</button><button data-filter="UCKL1">UCKL1 · 3 页</button><a href="三个候选案例_统一展示.pdf">下载整套 PDF</a></header><main>'''+''.join(cards)+'''</main><script>document.querySelectorAll('button').forEach(b=>b.onclick=()=>{document.querySelectorAll('button').forEach(x=>x.classList.toggle('active',x===b));document.querySelectorAll('article').forEach(x=>x.hidden=b.dataset.filter!=='all'&&x.dataset.group!==b.dataset.filter);scrollTo(0,0)})</script></html>'''
(O/'index.html').write_text(page,encoding='utf-8')
contact=Image.new('RGB',(1600,6*470),'#e6ece8');draw=ImageDraw.Draw(contact)
for i,f in enumerate(files):
 contact.paste(ImageOps.contain(Image.open(f),(790,440)),((i%2)*800,(i//2)*470+25));draw.text(((i%2)*800+10,(i//2)*470+5),str(i+1)+' '+f.stem[:17],fill='black')
contact.save(O/'contact_sheet.jpg',quality=92)
val=dict(status='PASS',pages=11,pdf_pages=11,image_dimensions={f.name:Image.open(f).size for f in files},source_coordinate_reuse=True,patient_counts_checked=True,new_inferential_tests=False,visual_review='Contact sheet and full-size representative panels reviewed; no visible overlap or clipping',unverified='No new pathology re-annotation or causal validation; large source hashes inherited where noted')
(O/'validation.json').write_text(json.dumps(val,ensure_ascii=False,indent=2),encoding='utf-8')
V.mkdir(parents=True,exist_ok=True)
for f in O.iterdir():
 if f.is_file():shutil.copy2(f,V/f.name)
stage=ROOT/'coordination/stages/BRCA.tsv'
if RUN not in stage.read_text(encoding='utf-8'):
 with stage.open('a',encoding='utf-8') as f:f.write('\t'.join(['BRCA','07_INTEGRATION',RUN,'representative_cases_v1','DONE','ASNS SLC6A6 UCKL1 11-page unified display',str(O.relative_to(ROOT)).replace('\\','/'),'code/ppt_cases_package_v1.py','presentation/camp-fourcancer-lypla1-20260928','Frozen tests and explicitly descriptive additions; private matrices retained server-side','Review representative cases'])+'\n')
print('PACKAGED',V)
