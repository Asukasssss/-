"""Publish a reader-facing catalog of existing single-cell results, no reanalysis."""
from pathlib import Path
import csv,json,hashlib,subprocess
R=Path(__file__).resolve().parents[1]
RUN='20260922T113500Z_sc_catalog_v1'
O=R/'results/BRCA/06_EXTERNAL'/RUN;O.mkdir(parents=True,exist_ok=True)
base=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()
url='https://github.com/Asukasssss/-/blob/'+base+'/'
groups=[
 ('原117基因：两肿瘤研究与正常背景','06_EXTERNAL/20260919T140105Z_scRNA117_v1'),
 ('原117基因：来源稳定性','04_ROBUSTNESS/20260919T151200Z_all117_robustness_v2'),
 ('新增39基因：表达来源、点图与热图','06_EXTERNAL/20260922T023000Z_sc39_v1'),
 ('全156基因：Wu分型背景','06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1'),
 ('全156基因：表达UMAP','06_EXTERNAL/20260922T082327Z_sc156_umap_v1'),
 ('全156基因：表达与细胞类型对照UMAP','06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1'),
 ('全156基因：上皮专用三联图','06_EXTERNAL/20260922T092000Z_epithelial_umap_v1'),
]
tracked=set(subprocess.check_output(['git','-c','core.quotepath=false','ls-files'],cwd=R,encoding='utf-8').splitlines())
rows=[];parts=['# BRCA单细胞结果：表格、图和逐基因UMAP总入口','','这里集中展示全部候选的既有单细胞表达来源结果。当前150个基因＋6个历史基因均保留；新目录只做链接与清单，不新增分析。所有链接固定到已上传提交，便于另一个账号直接使用。','','## 先打开这些图','']
pdfs=[('全156基因表达图','06_EXTERNAL/20260922T082327Z_sc156_umap_v1/all156_expression_UMAP.pdf'),('全156基因表达—细胞类型对照图','06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/all156_expression_celltype_pairs.pdf'),('全156基因上皮专用三联图','06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/all156_epithelial_triptychs.pdf'),('全156基因分型热图','06_EXTERNAL/20260922T064625Z_sc156_subtypes_v1/all156_subtype_expression_patterns.pdf')]
for name,p in pdfs:parts.append(f'- [{name}]({url}results/BRCA/{p})')
parts+=['','PDF有全部基因；重点基因另提供PNG。页码请查对应目录的gene_page_index.tsv或gene_scales_and_pages.tsv。GitHub若不能在线预览大PDF，可下载原文件查看。','','## 来源和图形的边界','','Wu沿用作者注释，Pal采用既有重注释；Reed为正常背景，不能计为第三套独立肿瘤验证。供者标签不是本次新核验的独立患者。分型结果以Wu为主，Pal未接入可靠分型连接。表达图显示RNA位置，不证明代谢物来源、功能或通量。上皮细胞状态不是临床患者亚型。本次不新增差异检验、聚类、通讯或拟时序。','']
for title,sub in groups:
 d=R/'results/BRCA'/sub
 parts+=['','## '+title,'',f'[说明与方法]({url}{d.relative_to(R).as_posix()}/README_CN.md)','','|文件|用途|','|---|---|']
 for f in sorted(d.iterdir()):
  if not f.is_file() or f.name.startswith('.'):continue
  if sub.startswith('04_') and ('composition' in f.name or f.name.startswith('all174')):continue
  p=f.relative_to(R).as_posix();assert p in tracked,p
  kind={'.pdf':'完整图册','.png':'图片','.tsv':'数值或来源表','.xlsx':'Excel汇总','.json':'参数或验证','.py':'复现脚本','.md':'中文说明'}.get(f.suffix,'校验/记录')
  rows.append(dict(section=title,path=p,kind=kind,bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),source_commit=base,url=url+p))
  parts.append(f'|[{f.name}]({url+p})|{kind}|')
parts+=['','## 重点基因图：直接查看','','这些只是方便浏览的示例，不是筛选后的唯一名单。全部156基因见上方PDF。','']
for gene in ['ASNS','LYPLA1','ABHD12','LYPLA2','ENPP2','LPCAT1','GLS','GPCPD1','GPI','NNMT','PYCR1','ACSL4']:
 p='results/BRCA/06_EXTERNAL/20260922T085000Z_sc156_umap_pairs_v1/'+gene+'_expression_celltype_pair.png'
 q='results/BRCA/06_EXTERNAL/20260922T092000Z_epithelial_umap_v1/'+gene+'_epithelial_triptych.png'
 assert p in tracked and q in tracked
 parts+=[f'- **{gene}**：[全细胞表达与类型对照]({url+p}) · [上皮三联图]({url+q})']
parts+=['','## 与统一候选表连接','','[BRCA全部结果入口](BRCA_CURRENT_RESULTS_CN.md)保留配对代谢物、患者相关、RNA及解释。不能因为只想看图片，就把缺测、相反结果或历史候选删掉。']
(R/'docs/BRCA_SINGLE_CELL_RESULTS_CN.md').write_bytes(('\n'.join(parts)+'\n').encode())
def out(name,obj):(O/name).write_bytes((json.dumps(obj,ensure_ascii=False,indent=2)+'\n').encode())
with (O/'single_cell_artifact_catalog.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
out('analysis_spec.json',dict(operation='catalog_only',source_commit=base,new_statistics=0,new_images=0,scope='All published BRCA single-cell source artifacts in seven registered runs; mixed composition-analysis files excluded'))
out('validation.json',dict(status='DONE',catalog_files=len(rows),groups=len(groups),all_files_tracked=True,all_paths_exist=True,all_files_sha256_recorded=True,full156_pdf_links=3,selected_gene_png_pairs=12,no_patient_matrix_read=True))
out('source_manifest.json',dict(source_commit=base,source_directories=[x[1] for x in groups],hashes='single_cell_artifact_catalog.tsv'))
(O/'README_CN.md').write_bytes(('''# BRCA单细胞结果集中交付

## 本轮问题
将单细胞数值、来源稳定性、分型、热图、点图与全156基因UMAP集中呈现。
## 输入与范围
7个已发布运行目录，源提交见analysis_spec。未读取患者矩阵。
## 实际结果
完整文件数见validation，逐文件哈希见single_cell_artifact_catalog.tsv。
## 新手解释
先看完整图册，再用页码索引找基因，最后对照来源表的供者覆盖与表达差距。
## 限制/反证
仅整理入口，无新统计或新图；表达来源不代表功能或代谢物来源。
## 当前决定
全156均保留，单细胞停在来源；重点PNG不是唯一候选。
## 下一步
从统一候选表讨论有限问题，不自动新增机制分析。
## 复现
python code/brca_single_cell_catalog_v1.py；再次发布用新RUN_ID。

[单细胞图表总入口](../../../../docs/BRCA_SINGLE_CELL_RESULTS_CN.md)
''').encode())
(O/'.gitattributes').write_bytes(b'* -text\n')
reg=R/'coordination/stages/BRCA.tsv'
with reg.open(encoding='utf-8') as f:a=list(csv.DictReader(f,delimiter='\t'))
if not any(x['run_id']==RUN for x in a):
 a.append(dict(cancer='BRCA',stage_id='06_EXTERNAL',run_id=RUN,analysis_version='sc_catalog_v1',status='DONE',scope=f'{len(rows)} single-cell files indexed;all156 PDF atlases;12 gene PNG pairs;7 source runs',result_path=O.relative_to(R).as_posix(),code_path='code/brca_single_cell_catalog_v1.py',git_branch='analysis/brca-functional-review-20260919',reason='Catalog only;no new images statistics or patient validation',next_action='Use visible single-cell entry alongside full candidate records'))
 a.sort(key=lambda x:(x['stage_id'],x['run_id']))
 with reg.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(a[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(a)
print(json.dumps(dict(files=len(rows),groups=len(groups))))
