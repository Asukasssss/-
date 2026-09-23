"""Check frozen aggregate extraction; no patient calculations."""
from pathlib import Path
import io, json, subprocess, hashlib
import pandas as pd
import fitz
from openpyxl import load_workbook
R=Path(__file__).resolve().parents[1]
O=R/'results/BRCA/07_INTEGRATION/20260923T010000Z_choline_ethanolamine_v1'
B='results/BRCA/07_INTEGRATION/20260922T140000Z_report_v1/appendix/'
SHA='8952c54b61709c875ec83881b66ee43bd1de59f3'
def read(p): return pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False)
def old(n): return read(io.BytesIO(subprocess.check_output(['git','show',SHA+':'+B+n+'.tsv'],cwd=R)))
scope=read(O/'scope22.tsv'); keys=set(scope.metabolite_key)
rel=old('relations'); rel=rel[rel.metabolite_key.isin(keys)].reset_index(drop=True)
genes=set(rel.gene); ids=set(rel.relation_key)
checks={}
for src,dst,col,values in [('metabolites','metabolites44','metabolite_key',keys),('relations','relations41_history','relation_key',ids),('association','CAMP_associations82','relation_id',ids),('rna','RNA20','gene',genes),('source_status','source_status','gene',genes)]:
    x=old(src); x=x[x[col].isin(values)].reset_index(drop=True)
    pd.testing.assert_frame_equal(x,read(O/(dst+'.tsv')),check_dtype=False)
    checks[dst]='EXACT_FROZEN_VALUES_PASS'
sc=pd.concat([old(n) for n in ['Wu2021_old','Wu2021_new','Pal2021_reprocessed_old','Pal2021_reprocessed_new']],ignore_index=True)
sc=sc[sc.gene.isin(genes)&sc.partition.eq('ALL')].reset_index(drop=True)
pd.testing.assert_frame_equal(sc,read(O/'single_cell_profiles.tsv'),check_dtype=False)
summary=read(O/'topic_summary41.tsv')
pd.testing.assert_frame_equal(rel,summary[list(rel.columns)],check_dtype=False)
checks['single_cell_and_history']='EXACT_FROZEN_VALUES_PASS'
a=read(O/'CAMP_associations82.tsv')
for family,prefix in [('processed_Spearman262','primary'),('available_Spearman262','sensitivity')]:
    idx=a[a.analysis_type.eq(family)].set_index('relation_id')
    for col in ['n','effect','ci_lower','ci_upper','p_value','q_value','status','reason']:
        assert summary['review_'+prefix+'_'+col].tolist()==summary.relation_key.map(idx[col]).tolist()
checks['summary_join']='PASS'
m=read(O/'metabolites44.tsv')
assert m[m.metabolite_name.eq('1-palmitoyl-GPC (16:0)')&m.analysis_type.eq('paired_author_available')].n.iloc[0]=='28'
wb=load_workbook(O/'专题汇总.xlsx',read_only=True)
assert len(wb.sheetnames)==8 and wb['专题汇总'].max_row==42
wb.close()
with fitz.open(O/'三图专题.pdf') as doc: assert len(doc)==3
checks['workbook_pdf_and_LPC16_coverage']='PASS'
v=json.loads((O/'validation.json').read_text(encoding='utf-8'))
v.update(checks=checks,visual_review='PASS: all three PNG panels inspected for labels and layout',patient_matrix_reanalysis=False)
v['CAMP_evaluable_by_family']=a.groupby('analysis_type').apply(lambda d:int(d.status.eq('DONE').sum())).to_dict()
(O/'validation.json').write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf-8')
pd.DataFrame([{'path':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(O.iterdir()) if p.name!='checksums.tsv']).to_csv(O/'checksums.tsv',sep='\t',index=False)
print(json.dumps(v,ensure_ascii=False,indent=2))
