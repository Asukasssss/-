"""Source admission diagnostics only; deliberately calculates no associations."""
import pathlib,csv,json,math,hashlib,collections,openpyxl,sys
root=pathlib.Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
run=pathlib.Path(sys.argv[1]).resolve();assert run.parent==root/'results/collaborative/COAD/B' and (run/'.running').is_dir()
old=root/'results/collaborative/COAD/B/20260921T050501Z_integral_acquire_v1'
src=root/'data/candidates/coad_integral_20260921/ac4c04421_si_003.xlsx'
assert hashlib.sha256(src.read_bytes()).hexdigest()=='30037314867bf80f4b05b3ced26766094aa6e63ae939bf45b69731819358d67c'
w=openpyxl.load_workbook(src,read_only=True,data_only=True)
meta=list(w['A_sample_information'].values)[1:];tumor=[x[0]for x in meta if x[2]=='Tumor'];assert len(tumor)==6
assert len({x[1]for x in meta if x[0]in tumor})==6
tabs={};diagnostics=[]
for sn,start,key in [('B_metabolome',1,0),('G_transcriptome_gene',2,0),('D_proteome_protein',2,1)]:
 rows=list(w[sn].values);h=rows[0];assert len(h[start:])==len(set(h[start:]))==12 and set(h[start:])=={x[0]for x in meta}
 ii=[h.index(t)for t in tumor];d=collections.defaultdict(list)
 for row in rows[1:]:d[str(row[key])].append([row[i]for i in ii])
 tabs[sn]=d
 cols=list(zip(*(row[start:]for row in rows[1:])));f=lambda v:isinstance(v,(float,int))and math.isfinite(v)
 sums=[sum(x for x in col if f(x))for col in cols];v=[x for col in cols for x in col if f(x)]
 diagnostics.append(dict(sheet=sn,records=len(rows)-1,exact_sample_label_set=True,column_order_matches_metabolome=h[start:]==list(w['B_metabolome'].values)[0][1:],finite_integer_fraction=sum(x==int(x)for x in v)/len(v),column_sum_min=min(sums),column_sum_max=max(sums),zero_cells=sum(x==0 for x in v),text_missing=sorted({str(x)for col in cols for x in col if not f(x)})))
with(old/'relation_coverage39.tsv').open()as f:relations=list(csv.DictReader(f,delimiter='\t'))
out=[]
for layer,sn in [('RNA','G_transcriptome_gene'),('PROTEIN','D_proteome_protein')]:
 for t in relations:
  gene=tabs[sn].get(t['gene'],[]);met=tabs['B_metabolome'].get(t['source_metabolite'],[])
  numeric=lambda v:len(v)==6 and all(isinstance(x,(int,float))and math.isfinite(x)for x in v)
  reason=[]
  if t['metabolite_match_status']!='NAME_COMPATIBLE':reason.append(t['metabolite_match_status'])
  if len(gene)!=1:reason.append('GENE_ROW_ABSENT_OR_NONUNIQUE')
  elif not numeric(gene[0]):reason.append('GENE_INCOMPLETE')
  elif len(set(gene[0]))==1:reason.append('GENE_CONSTANT')
  if len(met)!=1:reason.append('METABOLITE_ROW_ABSENT_OR_NONUNIQUE')
  elif not numeric(met[0]):reason.append('METABOLITE_INCOMPLETE')
  elif len(set(met[0]))==1:reason.append('METABOLITE_CONSTANT')
  out.append(dict(layer=layer,metabolite_key=t['metabolite_key'],metabolite_name=t['metabolite_name'],gene=t['gene'],source_metabolite=t['source_metabolite'],name_numeric_nonconstant_pass=not reason,structural_reason=';'.join(reason)or'PASS',gene_finite_n=sum(isinstance(v,(float,int))and math.isfinite(v)for v in gene[0])if len(gene)==1 else 'NA',gene_nonzero_n=sum(isinstance(v,(float,int))and math.isfinite(v)and v!=0 for v in gene[0])if len(gene)==1 else 'NA',analysis_status='NOT_RUN',admission_status='NEEDS_REVIEW',admission_reason='Cross-sample normalization and source missing/zero semantics not certified; no correlation computed'))
with(run/'admission78.tsv').open('w',newline='')as f:wr=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');wr.writeheader();wr.writerows(out)
summary=dict(status='PARTIAL',source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),patient_n=6,structural_pass={layer:sum(x['layer']==layer and x['name_numeric_nonconstant_pass']for x in out)for layer in ['RNA','PROTEIN']},diagnostics=diagnostics,source_admission='NEEDS_REVIEW',correlations_computed=0,zero_handling='Retained as numeric source zeros during diagnostics; semantics not certified',site_scope='CRC, anatomical site mix unverified')
(run/'admission_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
