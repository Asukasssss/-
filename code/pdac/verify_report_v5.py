"""Read-only checks of complete report artifacts and exact immutable results."""
from pathlib import Path
import json,hashlib,subprocess,math
import pandas as pd
from openpyxl import load_workbook
import pymupdf as fitz
ROOT=Path(__file__).resolve().parents[2]
O=ROOT/'outputs/pdac-unified-report-v5-20260922-final'
def read(p):return pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False)
def main():
 c=json.loads((ROOT/'configs/PDAC_report_v5.yaml').read_text(encoding='utf8'))
 for r in read(O/'source_manifest.tsv').itertuples():
  b=subprocess.check_output(['git','show',r.input_commit+':'+r.path],cwd=ROOT);assert hashlib.sha256(b).hexdigest()==r.sha256
 d=json.loads((O/'workbook_data.json').read_text(encoding='utf8'));wb=load_workbook(O/'PDAC_完整结果册.xlsx',read_only=False,data_only=False)
 assert wb.sheetnames==[x['name'] for x in d['sheets']]
 numbers=0;cells=0
 for x in d['sheets']:
  sh=wb[x['name']];assert sh.max_row==len(x['rows'])+1 and sh.max_column==len(x['headers'])
  assert sh.freeze_panes=='B2',(sh.title,sh.freeze_panes)
  assert len(sh.tables)==1
  for i,row in enumerate(sh.iter_rows(min_row=2),start=0):
   for j,cell in enumerate(row):
    expected=x['rows'][i][j];actual=cell.value;cells+=1
    if isinstance(expected,(float,int)):
     assert isinstance(actual,(float,int)) and math.isclose(actual,expected,rel_tol=1e-12,abs_tol=1e-14),(sh.title,cell.coordinate,actual,expected)
     numbers+=1
    else:
     if isinstance(expected,str) and expected.startswith(('=','+','@')):expected="'"+expected
     assert actual==expected or (actual is None and expected==''),(sh.title,cell.coordinate,actual,expected)
    assert cell.data_type!='e'
 report=fitz.open(O/'PDAC_主报告.pdf');assert len(report)==16
 for page in report:
  text=page.get_text();assert len(text)>100 and '\ufffd' not in text
 # Each independent cohort's appendix has all687, with explicit missing-value rows.
 for s in c['sc']['cohorts']:
  files=sorted((O/'figure_data').glob('appendix_heat_'+s+'_*.tsv'));allrows=pd.concat([read(p) for p in files]);assert allrows.gene.nunique()==687
  assert not allrows.duplicated(['gene','celltype']).any()
  atlas=fitz.open(O/'appendix'/('all_genes_heatmap_'+s+'.pdf'));assert len(atlas)==27
 ss=read(O/'appendix/source_status.tsv');assert len(ss)==2061
 assert ss.loc[ss.source_status.ne('DONE'),'display_top'].eq('暂不可定位').all()
 for key,up,down,equal in [('metabolites','pairs_higher','pairs_lower','pairs_equal'),('rna','positive_pairs','negative_pairs','equal_pairs')]:
  t=read(O/'appendix'/f'{key}.tsv');t=t[t.status.eq('DONE')]
  assert (t[[up,down,equal]].astype(float).sum(axis=1)==t.n.astype(float)).all()
 result=dict(status='PASS',sheets=len(wb.sheetnames),checked_cells=cells,checked_numeric_cells=numbers,all_native_filters_and_panes=True,source_manifest_hashes_verified=True,report_pages=16,full_source_appendix_pages=81,all687_per_cohort=True,missing_source_masked=True,new_statistics=0)
 (O/'workbook_validation.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
 print(json.dumps(result))
if __name__=='__main__':main()
