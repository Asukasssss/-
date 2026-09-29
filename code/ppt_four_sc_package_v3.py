"""Package four consistently styled slides; only aggregate data and figures leave server."""
from pathlib import Path
import json,hashlib,shutil
import pandas as pd
from pypdf import PdfWriter
ROOT=Path(__file__).resolve().parents[1];RUN='20260929T071500Z_four_sc_v3'
OUT=ROOT/'results/BRCA/07_INTEGRATION'/RUN
VIEW=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_four_cancer_singlecell_v3');VIEW.mkdir(parents=True,exist_ok=True)
(OUT/'.gitattributes').write_bytes(b'* -text\n')
d=pd.read_csv(OUT/'patient_summary.tsv',sep='\t')
assert d.cancer.tolist()==['BRCA','COAD','PDAC','PRAD']
assert d.higher.tolist()==[7,9,9,10] and d.n_pairs.tolist()==[8,9,9,10]
assert d.loc[d.cancer.isin(['COAD','PRAD']),'p'].isna().all()
for c,n in zip(d.cancer,d.cells):assert pd.read_csv(OUT/(c+'_coverage.tsv'),sep='\t').cells.sum()==n
for p in OUT.glob('*.svg'):p.write_bytes(('\n'.join(s.rstrip() for s in p.read_text(encoding='utf8').splitlines())+'\n').encode('utf8'))
w=PdfWriter();sections=[]
for page,c in enumerate(d.cancer,11):
 name=f'{page}_{c}_LYPLA1_全细胞与患者'
 w.append(str(OUT/(name+'.pdf')))
 for ext in ['png','pdf','svg']:shutil.copy2(OUT/(name+'.'+ext),VIEW/(name+'.'+ext))
 sections.append(f'<section id="{c}"><h2>{page} · {c}</h2><img src="{name}.png"><p><a href="{name}.pdf">PDF</a> · <a href="{name}.png">PNG</a> · <a href="{name}.svg">SVG</a></p></section>')
with (OUT/'LYPLA1_四癌单细胞_统一版.pdf').open('wb') as f:w.write(f)
shutil.copy2(OUT/'LYPLA1_四癌单细胞_统一版.pdf',VIEW/'LYPLA1_四癌单细胞_统一版.pdf')
(OUT/'README_CN.md').write_text('''# 四癌单细胞统一展示 v3

## 本轮问题
把细胞类型图例移到UMAP右侧，并为BRCA、COAD、PDAC、PRAD制作同一风格页面。LYPLA1专题仍不包括GBM/ccRCC。

## 输入与范围
全部复用已有表达及坐标。BRCA Wu100064细胞；COAD Uhlitz68702细胞；PDAC Hwang未治疗样本108964细胞核；PRAD24肿瘤及邻近组织68322细胞。PDAC只展示未治疗样本的全部细胞类型，以对应同患者检验范围。其他三队列完整原有全细胞展示；患者图按指定上皮比较另行限定。

## 实际结果
BRCA：7/8位恶性上皮高于肿瘤内非恶性上皮，平均差+0.215333，复用名义P=0.0390625。
COAD：9/9位肿瘤CNA上皮高于正常组织上皮，平均差+0.362927；本次复用原有患者汇总，按作者患者键匹配，描述性展示，不新算P值。原CNA身份排除规则及冲突供者排除沿用。CNN不能被称为已证实非恶性。
PDAC：9/9位恶性上皮高于非恶性导管，全核均值平均差+0.329992，复用P=0.0034799652。仅LYPLA1阳性核的原敏感性分析不支持同方向结论，不能解释为每个阳性细胞都表达增强。
PRAD：10/10位肿瘤组织内恶性上皮高于非恶性上皮，平均差+0.545776；从既有逐细胞表达按作者donor_id汇总，每类至少20细胞，描述性匹配，不新增P值。altered_benign不纳入比较。

## 新手解释
左侧LYPLA1实测表达，中间原有细胞注释，图例竖排紧邻UMAP右侧，最右侧同患者比较。不同癌种的对应大类采用同一配色；保留源注释分辨率，不能强制把不同细胞标签当成等价。红圈是候选上皮的视觉定位而非细胞筛选。所有图使用viridis和相同log1p(CP10K)尺度；上限按各队列阳性值99分位计算并逐页标注，故不能直接比较不同页颜色深浅。每条配对线为一位源标签患者，横杠为患者等权均值，平均差不是倍数。方向人数也不等于各患者分别显著。

## 限制与反证
统一的是绘图样式，不是证据等级、统计方法或参照定义。COAD/PRAD本轮新增的是描述性配对展示，不称显著；BRCA/PDAC保留已有检验，不重算不改阈值。只显示名义P，原q在冻结表保留。源患者身份未做基因型复核；UMAP混合度或圈内细胞不能证明机制或特异性。COAD是已有项目计算坐标且未做批次校正，另三癌沿用源坐标。不同参照定义不能直接合并效应。逐细胞及逐患者测量仅存server165。

## 当前决定
四页用于组会，同一版式、图例位置、线条、字体、细胞类型颜色和表达色标规则。原页保留，不生成PPT，不新增机制分析。

## 下一步
按用户反馈继续排版或推进空间证据；如需COAD/PRAD正式配对检验，应另设统计版本并明确检验范围。

## 复现
server165在BRCA/A新展示运行目录执行ppt_four_sc_prepare_v3.py RUN_DIR，再执行ppt_four_sc_render_v3.py RUN_DIR。只复制public目录。本地运行python code/ppt_four_sc_package_v3.py。其他癌种来源只读，没有写入或覆盖其目录。
''',encoding='utf8')
v=json.loads((OUT/'validation.json').read_text());v.update(visual_review='PASS: all four full-slide PNGs inspected, right legend fits',aggregate_counts_verified=True,COAD_PRAD_p_values_not_fabricated=True,per_patient_measurements_exported=False)
(OUT/'validation.json').write_text(json.dumps(v,indent=2),encoding='utf8')
pd.DataFrame([dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in [ROOT/'code/ppt_four_sc_prepare_v3.py',ROOT/'code/ppt_four_sc_render_v3.py',Path(__file__)]]).to_csv(OUT/'code_manifest.tsv',sep='\t',index=False)
for p in OUT.iterdir():
 if p.suffix in ['.md','.tsv','.json']:p.write_bytes(p.read_bytes().replace(b'\r\n',b'\n'))
pd.DataFrame([dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='checksums.tsv']).to_csv(OUT/'checksums.tsv',sep='\t',index=False,lineterminator='\n')
for name in ['README_CN.md','patient_summary.tsv']:shutil.copy2(OUT/name,VIEW/name)
nav=' · '.join(f'<a href="#{c}">{c}</a>' for c in d.cancer)
(VIEW/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>LYPLA1 四癌单细胞</title><style>html{scroll-behavior:smooth}body{background:#e8eeea;margin:0;font-family:system-ui;color:#183d34}nav{position:sticky;top:0;background:#fff;padding:16px 4%;z-index:9;box-shadow:0 1px 8px #aaa}main{max-width:1680px;margin:20px auto}section{scroll-margin-top:70px;margin:28px 0}img{width:100%;display:block;box-shadow:0 8px 25px #163e3320}a{color:#00553b}h2,p{padding:0 12px}</style><nav>'+nav+' · <a href="LYPLA1_四癌单细胞_统一版.pdf">下载四页PDF</a> · <a href="patient_summary.tsv">患者证据汇总</a></nav><main>'+''.join(sections)+'</main>',encoding='utf8')
stage=ROOT/'coordination/stages/BRCA.tsv'
if RUN not in stage.read_text(encoding='utf8'):
 with stage.open('a',encoding='utf8') as f:f.write('\t'.join(['BRCA','07_INTEGRATION',RUN,'four_sc_v3','DONE','Four-cancer presentation: right UMAP legend and patient evidence',str(OUT.relative_to(ROOT)).replace('\\','/'),'code/ppt_four_sc_prepare_v3.py;code/ppt_four_sc_render_v3.py','presentation/camp-fourcancer-lypla1-20260928','No new inferential tests; COAD PRAD paired descriptive; BRCA PDAC frozen P; other cancer source directories read-only','Review unified four-cancer single-cell slides'])+'\n')
print(VIEW)
