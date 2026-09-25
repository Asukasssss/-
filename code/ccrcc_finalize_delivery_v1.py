"""Package only ccRCC public summaries and code; verify hashes and repository bytes."""
import argparse,hashlib,json,subprocess,zipfile
from pathlib import Path
import pandas as pd
import openpyxl,pymupdf
R=Path(__file__).resolve().parents[1];C=R/'results/ccRCC'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def main(report_run):
 out=C/'07_INTEGRATION'/report_run
 rows=[]
 specs=[
 ('01_CAMP','20260925T143000Z_discovery_v1','NEEDS_REVIEW','ccRCC4 initial histology audit','superseded by original-study histology v2','ccrcc_discovery_v1.py'),
 ('02_MAPPING','20260925T143000Z_discovery_v1','NEEDS_REVIEW','initial ccRCC4 mapping','superseded workpool;keep historical genes','ccrcc_mapping_v1.py'),
 ('03_PATIENT','20260925T143000Z_discovery_v1','NEEDS_REVIEW','initial ccRCC4 patient results','superseded histology scope','ccrcc_patient_v1.py'),
 ('04_ROBUSTNESS','20260925T143000Z_discovery_v1','NEEDS_REVIEW','initial ccRCC4 paired discovery','superseded histology scope','ccrcc_discovery_v1.py'),
 ('01_CAMP','20260925T151000Z_ccrcc3_v1','DONE','ccRCC3 author identity;17 pairs','later therapy audit independently confirms histology and RNA links','ccrcc3_discovery_v1.py'),
 ('02_MAPPING','20260925T151000Z_ccrcc3_v1','PARTIAL','ccRCC3 221 relations;144 genes','401 workpool features NEEDS_REVIEW;bounded library','ccrcc_mapping_v1.py'),
 ('03_PATIENT','20260925T151000Z_ccrcc3_v1','DONE','ccRCC3 primary/available/RNA','coarse therapy sensitivity superseded separately','ccrcc3_patient_v1.py'),
 ('04_ROBUSTNESS','20260925T151000Z_ccrcc3_v1','DONE','ccRCC3 711 features;17pairs','mixed therapy;not independent replication','ccrcc3_discovery_v1.py'),
 ('01_CAMP','20260925T155000Z_ccrcc4_histology_v2','DONE','ccRCC4 corrected pathology;12pairs','5 specimens from3patients excluded;history preserved','ccrcc4_discovery_v2.py'),
 ('02_MAPPING','20260925T155000Z_ccrcc4_histology_v2','PARTIAL','ccRCC4 130 relations;95 current genes','388 workpool features NEEDS_REVIEW;bounded library','ccrcc_mapping_v1.py'),
 ('03_PATIENT','20260925T155000Z_ccrcc4_histology_v2','DONE','ccRCC4 corrected all-family statistics','12 patients;treated context','ccrcc4_patient_v2.py'),
 ('04_ROBUSTNESS','20260925T155000Z_ccrcc4_histology_v2','DONE','ccRCC4 corrected 904feature discovery','12 pairs;available sensitivity kept separately','ccrcc4_discovery_v2.py'),
 ('03_PATIENT','20260925T155500Z_ccrcc3_therapy_v2','DONE','ccRCC3 refined therapy;primary/RNA exact reuse','NO includes1cabozantinib;combined author strata','ccrcc3_treatment_sensitivity_v2.py'),
 ('05_FUNCTION',report_run,'NOT_RUN','no new functional/mechanistic analysis','outside current discovery-to-cell-source scope','ccrcc_integrate_report_v1.py'),
 ('06_EXTERNAL','20260925T153000Z_source_union_v1','PARTIAL','GSE159115 all157genes source;one study complete','second independent study/UMAP not completed','ccrcc_sc_source_v1.py'),
 ('07_INTEGRATION',report_run,'DONE','351relations;154current/157history genes;PDF Excel and ZIP','scientific gaps explicitly retained','ccrcc_integrate_report_v1.py')]
 branch=subprocess.check_output(['git','branch','--show-current'],cwd=R,text=True).strip()
 for stage,run,status,scope,reason,code in specs:
  result=C/stage/run if stage!='05_FUNCTION' else out
  rows.append(dict(cancer='ccRCC',stage_id=stage,run_id=run,analysis_version=code.removesuffix('.py'),status=status,scope=scope,result_path=result.relative_to(R).as_posix(),code_path='code/'+code,git_branch=branch,reason=reason,next_action='read current integration and explicit remaining gaps'))
 save(pd.DataFrame(rows).sort_values(['stage_id','run_id']),R/'coordination/stages/ccRCC.tsv')
 # Add only missing public artifacts to existing checksum inventories; do not change numerical tables.
 for checksum in C.rglob('checksums.tsv'):
  stage=checksum.parent
  if stage==out:continue
  files=[p for p in stage.rglob('*') if p.is_file() and p.name!='checksums.tsv']
  prior=pd.read_csv(checksum,sep='\t')
  for _,z in prior.iterrows():assert sha(stage/z['file'])==z.sha256,(stage,z['file'])
  save(pd.DataFrame([dict(file=p.relative_to(stage).as_posix(),sha256=sha(p)) for p in sorted(files)]),checksum)
 archive=out/'ccRCC_完整交付包.zip'
 doc=pymupdf.open(out/'ccRCC_分析报告.pdf');assert len(doc)==14
 text=[p.get_text() for p in doc];assert all(len(t)>70 and '\ufffd' not in t for t in text)
 wb=openpyxl.load_workbook(out/'ccRCC_完整分析与候选比较.xlsx',read_only=False,data_only=False)
 assert len(wb.sheetnames)==18
 for ws in wb:
  assert ws.freeze_panes=='A2' and ws.auto_filter.ref
  for row in ws:
   for cell in row:assert cell.data_type!='f','unexpected spreadsheet formula'
 wb.close()
 source=pd.read_csv(out/'source_manifest.tsv',sep='\t')
 for _,z in source.iterrows():
  p=R/z.path_or_url;assert sha(p)==z.sha256
  tracked=subprocess.check_output(['git','show',z.input_commit+':'+z.path_or_url],cwd=R)
  assert hashlib.sha256(tracked).hexdigest()==z.sha256,('Git input differs',z.path_or_url)
 independent=json.loads((out/'independent_numeric_validation.json').read_text());assert independent['status']=='DONE'
 validation=json.loads((out/'validation.json').read_text());validation.update(PDF_text_all14pages_checked=True,PDF_representative_pages_visually_checked=True,Excel_18_sheets_filter_freeze_checked=True,spreadsheet_no_formulas=True,all_report_input_hashes_match_pinned_Git=True,independent_all_feature_effect_counts_verified=True,patient_private_files_exported=False)
 (out/'validation.json').write_text(json.dumps(validation,indent=2),encoding='utf-8')
 (C/'README_CN.md').write_text(f'# ccRCC 当前交付入口\n\n当前报告：[14页中文PDF](07_INTEGRATION/{report_run}/ccRCC_分析报告.pdf)；[Excel](07_INTEGRATION/{report_run}/ccRCC_完整分析与候选比较.xlsx)；[完整ZIP](07_INTEGRATION/{report_run}/ccRCC_完整交付包.zip)。\n\n采用ccRCC3 17对、ccRCC4原研究病理修正版12对；351条关系、154当前/157历史并集基因。患者主关联没有条目通过各队列q<0.05。\n\n[版本修正与边界](VERSION_DECISIONS_CN.md)。映射尚有待审项；单细胞仅GSE159115一项研究，无独立来源复现、UMAP或功能机制验证。\n',encoding='utf-8')
 save(pd.DataFrame([dict(file=p.relative_to(out).as_posix(),sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file() and p.name not in ['checksums.tsv',archive.name]]),out/'checksums.tsv')
 # Explicit white list; no private/server/source matrices and no unrelated cancer results.
 files=[p for p in C.rglob('*') if p.is_file() and p!=archive]
 files += list((R/'code').glob('ccrcc*.py'))+list((R/'code').glob('ccrcc*.R'))
 files += [R/'coordination/stages/ccRCC.tsv',R/'docs/STAGE_STANDARD_CN.md',R/'docs/CAMP_DISCOVERY_TO_CELL_SOURCE_SOP_CN.md']
 assert not any(p.suffix.lower() in ['.h5','.h5ad','.rds','.rdata','.gz'] or 'private' in p.name.lower() for p in files)
 with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in sorted(set(files)):z.write(p,p.relative_to(R).as_posix())
 with zipfile.ZipFile(archive) as z:
  assert z.testzip() is None
  for p in files:assert hashlib.sha256(z.read(p.relative_to(R).as_posix())).hexdigest()==sha(p)
 # The ZIP embeds its pre-ZIP checksum list; its own checksum stays external to avoid self-reference.
 (out/'PACKAGE_SHA256.txt').write_text(sha(archive)+'  '+archive.name+'\n',encoding='utf-8')
 print(json.dumps(dict(status='DONE',zip_files=len(files),zip_bytes=archive.stat().st_size,zip_sha256=sha(archive),PDF_pages=14,Excel_sheets=18,input_commit=source.input_commit.iloc[0])))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--report-run',required=True);a=p.parse_args();main(a.report_run)
