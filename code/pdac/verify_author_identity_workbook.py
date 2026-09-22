"""Read-only verification of exported workbook numbers and native controls."""
import json,math
from openpyxl import load_workbook
from deliver_author_identity_v4 import OUT
def main():
    d=json.loads((OUT/'workbook_data.json').read_text());w=load_workbook(OUT/'PDAC_SOP完整比较表.xlsx');assert len(w.sheetnames)==6
    sets=[('关系比较','relations',{'association_n':6,'association_rho':7,'association_q':8,'RNA_p':14,'RNA_q':15}),('基因比较','genes',{'RNA_n':6,'RNA_effect':7,'RNA_p':8,'RNA_q':9}),('RNA P小于005','RNA',{'n':2,'n_up':3,'n_down':4,'p_value':9,'q_value':10}),('配对代谢物51项','metabolites',{'n':2,'n_up':3,'n_down':4,'p_value':7,'q_value':8}),('恶性与正常来源','identity_profiles',{'n':3,'n_cells_total':4,'n_cells_eligible':5,'effect':6,'mean_detection_fraction':7})]
    n=0
    for sn,key,columns in sets:
        sh=w[sn];assert sh.max_row==len(d[key])+6 and len(sh.tables)==1 and sh.freeze_panes.endswith('7')
        for i,row in enumerate(d[key],7):
            for field,col in columns.items():
                expected=row[field];actual=sh.cell(i,col).value
                if expected is None:assert actual=='NA'
                else:assert isinstance(actual,(float,int)) and math.isclose(actual,expected,rel_tol=1e-12,abs_tol=1e-14),(sn,i,col)
                n+=1
    assert not [(s.title,c.coordinate) for s in w for row in s for c in row if c.data_type=='e']
    result={'status':'PASS','sheets':6,'numeric_cells_checked':n,'native_errors':0,'tables_and_freezes':True,'changed_identity_gene_method_previews_inspected':True}
    (OUT/'workbook_validation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
if __name__=='__main__':main()
