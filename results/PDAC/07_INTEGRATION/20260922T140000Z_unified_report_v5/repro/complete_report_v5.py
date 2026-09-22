"""Complete public provenance, interpretation and acceptance for the PDAC report."""
from pathlib import Path
import argparse,json,io,hashlib,subprocess,shutil,os
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/pdac-unified-report-v5-20260922-final'
PUBLIC=ROOT/'results/PDAC/07_INTEGRATION/20260922T140000Z_unified_report_v5'
A='results/PDAC/07_INTEGRATION/20260922T135400Z_report_adapter_v5/'
S='results/PDAC/06_EXTERNAL/20260922T131800Z_author_identity_v4/'
def tsv(p,d):d.to_csv(p,sep='\t',index=False,lineterminator='\n',na_rep='NA')
def read(p):return pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False)
def delivery_files():
 for root,dirs,files in os.walk(OUT,followlinks=False):
  dirs[:]=[x for x in dirs if x not in {'node_modules','rendered_pages','workbook_previews'}]
  for name in files:
   if name not in {'workbook_data.json','build.mjs','checksums.tsv','package_validation.json'} and not name.endswith('.inspect.ndjson'):yield Path(root)/name
def extras():
 c=json.loads((ROOT/'configs/PDAC_report_v5.yaml').read_text(encoding='utf8'));manifest=read(OUT/'source_manifest.tsv')
 more={
 'frozen_CAMP303':'results/PDAC/01_CAMP/20260920T051000Z_frozen_v1/results.tsv',
 'frozen_CAMP42':'results/PDAC/01_CAMP/20260920T051000Z_frozen_v1/significant_features.tsv',
 'direct_mapping':'results/PDAC/02_MAPPING/20260922T112600Z_sop_v3_mapping/direct_relations.tsv',
 'conditional_mapping':'results/PDAC/02_MAPPING/20260922T112600Z_sop_v3_mapping/conditional_relations.tsv',
 'mapping_evidence':'results/PDAC/02_MAPPING/20260922T112600Z_sop_v3_mapping/mapping_evidence.tsv',
 'gene_membership':'results/PDAC/02_MAPPING/20260922T112600Z_sop_v3_mapping/gene_pool_history_union.tsv',
 'metabolite_direction_ranking':'results/PDAC/01_CAMP/20260922T121000Z_sop_v3_discovery_complete/metabolite_direction_ranking.tsv',
 'author_identity_evidence':S+'author_identity_evidence.tsv',
 'malignant_normal_profiles':'results/PDAC/07_INTEGRATION/20260922T131900Z_author_identity_v4/author_malignant_normal_gene_profiles.tsv',
 'historical_relations':'results/PDAC/07_INTEGRATION/20260922T131900Z_author_identity_v4/historical_relation_index.tsv',
 'design_summary':A+'design_summary.tsv'}
 extra=[];new=[]
 for key,path in more.items():
  b=subprocess.check_output(['git','show',c['input_commit']+':'+path],cwd=ROOT)
  (OUT/'appendix'/f'{key}.tsv').write_bytes(b)
  extra.append(dict(key=key,path=path,input_commit=c['input_commit'],sha256=hashlib.sha256(b).hexdigest(),rows=len(pd.read_csv(io.BytesIO(b),sep='\t'))))
  if key=='direct_mapping' or key=='conditional_mapping':new.append(pd.read_csv(io.BytesIO(b),sep='\t',dtype=str,keep_default_na=False))
 manifest=pd.concat([manifest,pd.DataFrame(extra)],ignore_index=True).drop_duplicates('path');tsv(OUT/'source_manifest.tsv',manifest)
 wb=json.loads((OUT/'workbook_data.json').read_text(encoding='utf8'))
 profiles=pd.concat([read(OUT/'appendix'/f'{cohort}.tsv') for cohort in ['GSE263733','GSE278688','GSE242230']],ignore_index=True)
 first=['gene','cohort','celltype','effect','mean_detection_fraction','n','status']
 profiles=profiles[first+[col for col in profiles if col not in first]]
 for sheet in wb['sheets']:
  if sheet['name']=='单细胞来源全量':sheet.update(headers=list(profiles),rows=profiles.values.tolist(),numeric=[])
 mp=pd.concat(new,ignore_index=True)
 if not any(s['name']=='直接与条件映射' for s in wb['sheets']):wb['sheets'].append(dict(name='直接与条件映射',headers=list(mp),rows=mp.values.tolist(),numeric=[]))
 for sheet in wb['sheets']:
  for j,col in enumerate(sheet['headers']):
   if col.endswith('_id') or col.endswith('_key') or col in ['gene','input_commit','sha256','relation_ids']:continue
   vals=[row[j] for row in sheet['rows']];nonmissing=[v for v in vals if v not in ['NA','',None]]
   if not nonmissing:continue
   nums=pd.to_numeric(pd.Series(nonmissing),errors='coerce')
   if nums.notna().all():
    if col not in sheet['numeric']:sheet['numeric'].append(col)
    for row in sheet['rows']:
     if row[j] not in ['NA','',None]:row[j]=float(row[j])
 (OUT/'workbook_data.json').write_text(json.dumps(wb,ensure_ascii=False),encoding='utf8')
 rel=read(OUT/'appendix/all_relations.tsv');rows=[]
 special={x['relation_id']:x for x in c['examples']}
 for r in rel.itertuples():
  e=special.get(r.relation_id)
  rows.append(dict(relation_id=r.relation_id,metabolite=r.metabolite_name,gene=r.gene,discussion_depth='SPECIFIC_DISPLAY_DISCUSSION' if e else 'GENERAL_RETENTION_NOTE',observed=f'metabolite P={r.metabolite_p};q={r.metabolite_q};association rho={r.association_rho};P={r.association_p};q={r.association_q};RNA P={r.RNA_p};q={r.RNA_q}',inference=e['reason'] if e else '保留该精确生化关系作为探索候选；未针对本行开展深入机制解释',counterevidence_or_gap=(e['limitation']+'；' if e else '')+r.data_limitations,cannot_infer='RNA丰度、相关和首位表达均不能证明酶活/通量、因果方向或唯一作用细胞',minimal_question='该精确关系的患者效应与来源是否在明确独立单位下稳定？本轮不新增机制实验',source='appendix/all_relations.tsv'))
 tsv(OUT/'interpretation_by_relation.tsv',pd.DataFrame(rows))
 lines=['# PDAC 事实、解释、反证与缺项','',f'输入：{c["input_commit"]}。全部715关系保留；5个示例为专门讨论，其余明确为通用保留说明。','', '## 直接观察到的数值','11个作者配对：307代谢物中51个P<0.05、0个q<0.05。21个不同作者病例键的直接关联354/357可评估，65个P<0.05、1个q<0.05。当前RNA247/250可评估、40个q<0.05。','', '## 可提出的解释','ABAT与GABA在内部主/可用值均有关联；AASS的患者支持弱；TDO2有稳定CAF首位来源。MGLL和SLC6A19在本轮扩大检验族后的q分别0.059与0.0531，不得借旧q声称本轮过FDR。','', '## 不能推断的内容','51项入口均未过代谢物FDR。来源不构成关系外部验证；最高表达不是唯一作用细胞；两套导管标签不等于恶性。三队列严格稳定来源72基因不是72个机制或治疗靶点。','', '## 下一步最小问题','先按精确关系ID阅读主/可用值和RNA各层结果。需要恶性来源结论时先解决另外两队列逐细胞作者标签；本报告不追加机制分析。','']
 for e in c['examples']:lines+=['### '+e['gene']+' / '+e['relation_id'],e['reason']+'。限制：'+e['limitation']+'。','']
 (OUT/'INTERPRETATION_CN.md').write_text('\n'.join(lines),encoding='utf8')
 for filename in ['progress_summary.tsv','module_gaps.tsv','sample_audit_summary.tsv','design_summary.tsv']:
  shutil.copyfile(ROOT/A/filename,OUT/filename)
 for p in ['code/pdac/prepare_report_v5.py','code/pdac/configure_report_v5.py','code/pdac/render_report_v5.py','code/pdac/build_report_v5_workbook.mjs','code/pdac/complete_report_v5.py']:
  shutil.copyfile(ROOT/p,OUT/'repro'/Path(p).name)
 (OUT/'repro/README_CN.md').write_text('在仓库中保留固定输入提交，安装requirements.txt并提供中文字体。只读：python code/pdac/render_report_v5.py --config configs/PDAC_report_v5.yaml --validate-only。生成：同命令加--out <全新目录>。随后complete_report_v5.py extras补完整公开附录；build_report_v5_workbook.mjs使用Codex bundled @oai/artifact-tool生成工作簿。最后complete_report_v5.py finalize校验打包。后处理OUT/PUBLIC常量须指向新的目录。\n',encoding='utf8')
 print('EXTRAS_READY')

def finalize():
 import sys
 sys.path.insert(0,str(ROOT/'code'));import camp_results_report as rep
 assert (OUT/'PDAC_完整结果册.xlsx').exists()
 v=json.loads((OUT/'validation.json').read_text(encoding='utf8'))
 assert json.loads((OUT/'visual_review.json').read_text(encoding='utf8'))['status']=='PASS'
 assert json.loads((OUT/'workbook_validation.json').read_text(encoding='utf8'))['status']=='PASS'
 v.update(status='DONE_PRESENTATION_WITH_DOCUMENTED_RESEARCH_GAPS',visual_review='PASS',workbook_native_validation='PASS',new_statistics=0,clinical_identity='NEEDS_REVIEW',malignant_identity='PARTIAL',subtypes='NOT_RUN',UMAP='NOT_RUN')
 rep.dump(OUT/'validation.json',v)
 for filename in ['build_report_v5_workbook.mjs','complete_report_v5.py','verify_report_v5.py','refine_report_v5_layout.py']:
  shutil.copyfile(ROOT/'code/pdac'/filename,OUT/'repro'/filename)
 (OUT/'README_CN.md').write_text('''# PDAC CAMP至单细胞来源统一报告 v5

## 本轮问题
按交接规范接续已有结果，补真实报告与完整交付包；未重跑患者矩阵，未新增P/q。
## 输入与范围
输入提交2971100b37e40d66a88cf89102c0fd3a56b7bcee；统计基线bcccc21ff60d6c57111efe72804d50db8020245e；交接规范45c6ae9ab489a9bd4ead4e47d28185785a448156。
27肿瘤/12正常标本；11个作者明确配对用于代谢物和RNA；21个不同作者病例键用于肿瘤内关联。另6个肿瘤未默认独立。
## 实际结果
307项配对代谢物：51项P<0.05、0项q<0.05；原有值敏感性193项可评估、27项P<0.05、0项q<0.05。
357直接关系：354主分析可评估、65项P<0.05、1项q<0.05。358条件关系另族350项可评估、4项q<0.05。当前250基因RNA中247可评估、91项P<0.05、40项q<0.05。
全部687基因在三套单细胞均保留；可解释来源621/593/610个。三队列同首位且保持率均≥80%为72基因，当前直接池30个。
16页主PDF、Markdown、16页签完整Excel；97张新图各有PNG/PDF/SVG/图源/图注；3套各27页全基因热图册；15页既有全基因点图原PNG/SVG复用；完整ZIP可解压及逐文件哈希核对。
## 新手解释
先读PDAC_主报告.pdf，再按relation_id查appendix/all_relations.tsv或Excel完整关系页。357直接主关系与358条件补充不可合成一个q族。51入口的去向是30直接、12条件、3身份待复核、6本次无可靠精确关系；其中2个身份待复核特征仍有条件关系，所以实际没有计划关系的特征为7个。
## 限制/反证
代谢物51项均未过FDR。AASS患者统计未支持；MGLL与SLC6A19不能挪用早期候选池q。GSE242230拆分作者恶性和正常上皮；另两研究导管身份未定。临床身份/跨研究去重未认证，临床分型及UMAP未接入。源图没有逐患者或逐细胞原始值。
## 当前决定
本轮来源与报告范围已实际交付，研究缺项仍按PARTIAL/NOT_RUN/NEEDS_REVIEW记录；不声称整癌种研究完成或同关系外部验证。
## 下一步
按完整证据表讨论具体候选。若要恶性来源结论，需先取得另外两队列可连接的作者细胞判定；本轮不追加机制分析。
## 复现
python code/pdac/render_report_v5.py --config configs/PDAC_report_v5.yaml --validate-only
python code/pdac/render_report_v5.py --config configs/PDAC_report_v5.yaml --out <全新目录>
后处理与Artifact Tool工作簿步骤见repro/README_CN.md。公开PDF/Markdown/PNG/SVG及TSV入GitHub；XLSX/ZIP按仓库规定在本地交付。PR #2仍草稿、未合并。
''',encoding='utf8')
 acc=read(ROOT/'docs/PDAC/report_v5_reference/cancer_to_source_acceptance.tsv')
 for i,r in acc.iterrows():
  state='DONE';why='Actual output and immutable aggregate checks completed';path='validation.json'
  if r.check_id=='01':state='PARTIAL';why='Author pairs/keys and cross-omics matched;clinical identity/aliquot not recertified';path='sample_audit_summary.tsv'
  elif r.check_id=='07':state='PARTIAL';why='All687 covered;two cohorts lack exact-cell author malignant labels';path='appendix/author_identity_evidence.tsv'
  elif r.check_id=='09':state='PARTIAL';why='All/shared category summaries provided;clinical cross-study independence unresolved';path='appendix/cross_study.tsv'
  elif r.check_id=='10':state='NOT_RUN';why='No reliable clinical subtype join used;absence documented';path='module_gaps.tsv'
  elif r.check_id=='14':state='PARTIAL';why='Public package ready;remote publication must be verified after commit';path='PUBLICATION_CN.md'
  acc.loc[i,['status','evidence_path','reason']]=[state,path,why]
 tsv(OUT/'acceptance.tsv',acc)
 (OUT/'PUBLICATION_CN.md').write_text('本包固定输入2971100b37e40d66a88cf89102c0fd3a56b7bcee，数据基线bcccc21ff60d6c57111efe72804d50db8020245e。发布分支analysis/pdac-initial，PR #2草稿未合并。最终交付提交以Git历史及对话中的远程核验记录为准，避免自引用提交号。Excel和ZIP按仓库规则只保留本地交付；GitHub保留PDF、Markdown、图源及代码，原始患者/细胞数据不在包内。\n',encoding='utf8')
 # Explicit complete delivery; large machine tables split only in public copy.
 PUBLIC.mkdir(parents=True,exist_ok=False)
 excluded={'workbook_data.json'}
 for p in delivery_files():
  rel=p.relative_to(OUT)
  if p.name in excluded or rel.parts[0] in {'rendered_pages','workbook_previews'} or p.suffix in {'.xlsx','.zip'}:continue
  if rel.parts[0]=='figures' and p.suffix=='.pdf':continue # PNG+SVG and complete local PDF set retained
  if p.stat().st_size>5_000_000:
   if p.suffix=='.tsv':
    df=read(p)
    for i in range(0,len(df),3000):
     dst=PUBLIC/rel.parent/(p.stem+f'_part{i//3000+1:02d}.tsv');dst.parent.mkdir(parents=True,exist_ok=True);tsv(dst,df.iloc[i:i+3000])
    continue
   raise ValueError('Oversized public artifact '+str(rel))
  dest=PUBLIC/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
 # The machine JSON is an authoring intermediate, not a delivered source.
 files=list(delivery_files())
 import zipfile
 checks=pd.DataFrame([dict(path=p.relative_to(OUT).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(files)])
 tsv(OUT/'checksums.tsv',checks)
 zip_path=OUT.parent/'PDAC_统一报告完整包_v5.zip'
 with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as archive:
  for p in files+[OUT/'checksums.tsv']:archive.write(p,'PDAC_v5/'+p.relative_to(OUT).as_posix())
 with zipfile.ZipFile(zip_path) as archive:
  assert archive.testzip() is None
  for r in checks.itertuples():assert hashlib.sha256(archive.read('PDAC_v5/'+r.path)).hexdigest()==r.sha256
 rep.dump(OUT/'package_validation.json',dict(status='PASS',files=len(checks)+1,zip_sha256=hashlib.sha256(zip_path.read_bytes()).hexdigest(),zip_bytes=zip_path.stat().st_size))
 shutil.copyfile(OUT/'package_validation.json',PUBLIC/'package_validation.json')
 pubfiles=[p for p in PUBLIC.rglob('*') if p.is_file() and p.name!='checksums.tsv']
 tsv(PUBLIC/'checksums.tsv',pd.DataFrame([dict(path=p.relative_to(PUBLIC).as_posix(),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(pubfiles)]))
 print(json.dumps(dict(public=str(PUBLIC),zip=str(zip_path),files=len(files)),ensure_ascii=False))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['extras','finalize']);a=ap.parse_args();globals()[a.mode]()
