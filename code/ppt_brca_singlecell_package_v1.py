"""Package public slide artifacts only, after server rendering and visual inspection."""
from pathlib import Path
import shutil,json,hashlib
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
RUN='20260929T064100Z_brca_sc_slide_v1'
OUT=ROOT/'results/BRCA/07_INTEGRATION'/RUN
VIEW=Path('D:/CodexData/visualizations/2026/09/25/01a0d8cd-4545-75b0-ba43-922588deea1a/CAMP_BRCA_singlecell_v1');VIEW.mkdir(parents=True,exist_ok=True)
name='11_BRCA_LYPLA1_单细胞'
svg=OUT/(name+'.svg')
svg.write_bytes(('\n'.join(x.rstrip() for x in svg.read_text(encoding='utf8').splitlines())+'\n').encode('utf8'))
(OUT/'.gitattributes').write_bytes(b'* -text\n')
old=ROOT/'results/BRCA/06_EXTERNAL/20260923T040012Z_epithelial_paired20_v1/results.tsv'
d=pd.read_csv(OUT/'results.tsv',sep='\t'); o=pd.read_csv(old,sep='\t');o=o[(o.cohort=='Wu2021')&(o.gene=='LYPLA1')]
pd.testing.assert_frame_equal(d.reset_index(drop=True),o.reset_index(drop=True),check_dtype=False)
(OUT/'README_CN.md').write_text('''# BRCA 单细胞结果展示（第11页）

## 本轮问题
从CAMP组织表达转向细胞定位：LYPLA1是否在恶性上皮中高于同患者非恶性上皮？

## 输入与范围
Wu2021 / GSE176078。左、中图使用既有上皮嵌入全部28,844细胞，恶性24,489、非恶性4,355；坐标来自20260922T092000Z_epithelial_umap_v1，原有细胞身份保留。不是全细胞类型图，不声称恶性上皮在所有细胞类型中表达最高。右图引用20260923T040012Z_epithelial_paired20_v1既有配对结果，每患者每类至少20细胞，8位患者。UMAP展示范围大于配对检验范围。

## 实际结果
7/8位患者恶性上皮均值更高；平均患者内差0.2153329622，双侧精确配对符号秩P=0.0390625。作者Naive子集6/7同向，P=0.078125。图中标注名义P，按用户要求不以FDR作为展示门槛；原q完整保存在results.tsv。

## 新手解释
A说明哪些细胞被标注为恶性/非恶性；B在相同坐标显示LYPLA1实测RNA；C每条线对应一位患者，两端为该患者两类细胞的平均log1p(CP10K)，深绿短横为患者均值的组均值。平均差不是倍数。每个方向人数也不是单独通过显著检验的患者人数。

## 限制与反证
参照为肿瘤标本内非恶性上皮，不是健康人乳腺；患者身份沿用作者标签。UMAP复用既有分析者生成上皮嵌入，未做批次校正；未新聚类/重新判定恶性。图形不用于判定统计显著。单一队列、8位患者、后选择探索；作者Naive子集未达到名义显著。Pal现有子集缺非恶性上皮，不能作为该比较的复现。表达上调不等于酶活增强。

## 当前决定
用统一SYSU版式展示BRCA现有证据，不生成PPT。主结论限定为该队列恶性上皮上调的支持。

## 下一步
由用户选择其他癌种同类证据或BRCA空间定位；不自动增加新分析。

## 复现命令
在server165新运行目录复制字体、真实校标和冻结results.tsv，运行python3 ppt_brca_singlecell_slide_v1.py RUN_DIR。只回传public目录；逐细胞矩阵、坐标、患者表不回传。本地运行python code/ppt_brca_singlecell_package_v1.py完成打包。所有数值复用，无新检验。
''',encoding='utf8')
v=json.loads((OUT/'validation.json').read_text());v.update(local_frozen_results_identical=True,visual_review='PASS: full slide inspected; no clipping or overlap',public_only=True)
(OUT/'validation.json').write_text(json.dumps(v,indent=2),encoding='utf8')
for ext in ['png','pdf','svg']:shutil.copy2(OUT/(name+'.'+ext),VIEW/(name+'.'+ext))
for n in ['README_CN.md','results.tsv']:shutil.copy2(OUT/n,VIEW/n)
(VIEW/'index.html').write_text(f'<!doctype html><meta charset="utf-8"><title>CAMP · BRCA单细胞</title><style>body{{background:#e8eeea;margin:0;font-family:system-ui}}main{{max-width:1440px;margin:24px auto}}img{{width:100%;display:block;box-shadow:0 8px 25px #163e3320}}a{{color:#00553b}}</style><main><h2>第11页 · BRCA单细胞证据</h2><img src="{name}.png"><p><a href="{name}.pdf">高清PDF</a> · <a href="results.tsv">已有统计</a> · <a href="README_CN.md">方法与解释</a></p></main>',encoding='utf8')
pd.DataFrame([dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in [ROOT/'code/ppt_brca_singlecell_slide_v1.py',Path(__file__),old]]).to_csv(OUT/'local_source_manifest.tsv',sep='\t',index=False)
pd.DataFrame([dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='checksums.tsv']).to_csv(OUT/'checksums.tsv',sep='\t',index=False)
stage=ROOT/'coordination/stages/BRCA.tsv'
if RUN not in stage.read_text(encoding='utf8'):
 with stage.open('a',encoding='utf8') as f:f.write('\t'.join(['BRCA','07_INTEGRATION',RUN,'brca_sc_slide_v1','DONE','BRCA epithelial UMAP and frozen paired donor comparison',str(OUT.relative_to(ROOT)).replace('\\','/'),'code/ppt_brca_singlecell_slide_v1.py','presentation/camp-fourcancer-lypla1-20260928','No new tests; nonmalignant epithelial reference within tumor; 7 of 8 higher, nominal P=0.0391','Review BRCA single-cell slide'])+'\n')
print(VIEW)
