"""Structural read only. Do not print individual measurement rows."""
import pathlib,json,subprocess,openpyxl
root=pathlib.Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
dest=root/'data/candidates/coad_integral_20260921'
run=root/'results/collaborative/COAD/B/20260921T050501Z_integral_acquire_v1'
assert (run/'source_files.json').exists()
w=openpyxl.load_workbook(dest/'ac4c04421_si_003.xlsx',read_only=True,data_only=True)
summary=[]
for s in w:
 heads=list(s.iter_rows(max_row=2,values_only=True))
 summary.append(dict(sheet=s.title,rows=s.max_row,columns=s.max_column,headers=[[str(x)if x is not None else '' for x in row]for row in heads]))
(run/'structure_internal.json').write_text(json.dumps(summary,indent=2))
print(json.dumps([dict(sheet=s['sheet'],rows=s['rows'],columns=s['columns'])for s in summary]))
subprocess.run(['pdftotext','-layout',str(dest/'ac4c04421_si_001.pdf'),str(dest/'methods.txt')],check=True)
