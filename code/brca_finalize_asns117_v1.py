"""Append aggregate ASNS response evidence without modifying historical columns."""
from pathlib import Path
import csv,json,hashlib
ROOT=Path(__file__).resolve().parents[1]
RUN='20260921T034354Z_asns117_programs_v1'
OUT=ROOT/'results/BRCA/05_FUNCTION'/RUN
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:
  r=csv.DictReader(f,delimiter='\t');return r.fieldnames,list(r)
def write(p,c,r):
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,c,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(r)
def main():
 prior=ROOT/'results/BRCA/06_EXTERNAL/20260920T090922Z_all174_external117_resources_v1/all117_external_resources_appended.tsv'
 cols,rows=read(prior);assert len(rows)==117
 _,mapping=read(OUT/'orthology117.tsv');m={r['human_gene']:r for r in mapping}
 _,summary=read(OUT/'all117_two_construct_summary.tsv');s={(r['human_gene'],r['site']):r for r in summary}
 added=['Asns_response_mouse_gene','Asns_response_mapping_status','Asns_response_mapping_reason','Asns_response_Tumor','Asns_response_Lung','Asns_response_scope','Asns_response_source']
 assert not set(added)&set(cols)
 for r in rows:
  g=r['gene'];r.update(dict(zip(added,[m[g]['mouse_gene'],m[g]['status'],m[g]['reason'],s[g,'Tumor']['response_state'],s[g,'Lung']['response_state'],'mouse_selected_recultured_RNA;response_to_Asns_not_self_gene_function',OUT.relative_to(ROOT).as_posix()])))
 write(OUT/'all117_comparison_asns_response_appended.tsv',cols+added,rows)
 _,original=read(prior);assert [{c:r[c] for c in cols} for r in rows]==original
 (OUT/'integration_validation.json').write_bytes((json.dumps(dict(all117_retained=True,historical_columns_unchanged=len(cols),new_columns=len(added),prior_sha256=hashlib.sha256(prior.read_bytes()).hexdigest()),indent=2)+'\n').encode())
 sf,sr=read(ROOT/'coordination/stages/BRCA.tsv');sr=[r for r in sr if r['run_id']!=RUN]
 for stage,status,scope in [('05_FUNCTION','PARTIAL','117 retained;111 one2one orthologs;6 unresolved;3 fixed programs BH12;response not self-gene validation'),('07_INTEGRATION','DONE','All117 history retained;22 non-target Tumor concordant responders;no non-target Lung concordant responder')]:
  sr.append(dict(zip(sf,['BRCA',stage,RUN,'asns117_programs_v1',status,scope,OUT.relative_to(ROOT).as_posix(),'code/brca_finalize_asns117_v1.py','analysis/brca-functional-review-20260919','Mouse selected recultured cells;shared controls;no interaction test;fixed correlation0.01','Freeze this batch;compare appropriate direct GLS intervention;do not rerank all117 by response alone'])))
 write(ROOT/'coordination/stages/BRCA.tsv',sf,sorted(sr,key=lambda r:(r['stage_id'],r['run_id'])))
 manifests=[]
 for name in ['orthology_source_manifest.tsv','program_source_manifest.tsv']:
  _,rr=read(OUT/name);manifests.extend(rr)
 write(OUT/'source_manifest.tsv',['path','sha256'],manifests)
 scripts=['brca_asns117_map_v1.py','brca_asns_programs_v1.R','brca_asns117_validate_plot_v1.py','brca_finalize_asns117_v1.py']
 write(OUT/'code_manifest.tsv',['path','sha256'],[dict(path='code/'+s,sha256=hashlib.sha256((ROOT/'code'/s).read_bytes()).hexdigest()) for s in scripts])
 for p in OUT.iterdir():
  if p.suffix in {'.md','.tsv','.json','.txt'}:p.write_bytes(p.read_bytes().replace(b'\r\n',b'\n'))
 ch=OUT/'checksums.sha256';ch.write_bytes(''.join(f'{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n' for p in sorted(OUT.iterdir()) if p.is_file() and p!=ch).encode())
 print((OUT/'integration_validation.json').read_text())
if __name__=='__main__':main()
