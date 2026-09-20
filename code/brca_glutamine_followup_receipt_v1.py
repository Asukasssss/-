"""Reuse aggregate statistics and inspect public resource metadata; no patient reanalysis."""
from pathlib import Path
import csv,json,hashlib,requests,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
RUN='20260920T104956Z_glutamine_context_v1'
OUT=ROOT/'results/BRCA/04_ROBUSTNESS'/RUN
OLD=ROOT/'results/BRCA/06_EXTERNAL/20260920T090922Z_all174_external117_resources_v1'
def read(p):
 with p.open(encoding='utf-8-sig',newline='') as f:
  r=csv.DictReader(f,delimiter='\t');return r.fieldnames,list(r)
def write(p,fields,rows):
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 cols,rows=read(OLD/'external_all174_with_CAMP_context.tsv')
 primary=[r for r in rows if r['test_family']=='FUSCC_PRIMARY174']
 sig=[r for r in primary if r['status']=='DONE' and float(r['q_value'])<.05]
 opposite=[r for r in sig if r['direction_vs_CAMP_unadjusted']=='OPPOSITE_OR_ZERO']
 assert len(opposite)==8 and all(float(r['CAMP_unadjusted_q'])>=.05 for r in opposite)
 same=[r for r in sig if r['direction_vs_CAMP_unadjusted']=='SAME_SIGN']
 dual=[r for r in same if float(r['CAMP_unadjusted_q'])<.05]
 assert len(dual)==2 and {r['gene'] for r in dual}=={'ASNS','NNMT'}
 selected=[r for r in primary if (r['metabolite_key']=='KEGG:C00064' and r['gene'] in ['ASNS','GLS','GLS2','GLUL']) or r['gene']=='ASNS' or r in sig]
 write(OUT/'existing_relation_context.tsv',cols,selected)
 genes=['ASNS','GLS','GLS2','GLUL']
 models=[{'gene':g,'covariates':[h for h in genes if h!=g],'model':'other_three'} for g in genes]+[{'gene':'ASNS','covariates':['GLS'],'model':'mutual_pair'},{'gene':'GLS','covariates':['ASNS'],'model':'mutual_pair'}]
 spec={'analysis_version':'glutamine_context_v1','status':'NOT_RUN','planned_before_new_conditional_results':True,'selection':'post_selection_exploratory_using_existing_FUSCC_results;not_confirmatory','metabolite_key':'KEGG:C00064','external_peak':'M147T365_POS','genes':genes,'models':models,'unit':'author_patient_id','same_subset':'complete cases across all four genes and metabolite; no new imputation','method':'rank each variable; OLS residualize gene and metabolite against intercept and ranked covariates; Pearson correlation of residuals','p_method':'two-sided partial-correlation t approximation; df=n-k-2; not exact rank test','multiple_testing':'BH across six planned primary hypotheses; missing P internal1 and publicNA; Holm6 additionally reported','bootstrap':'2000 patient-row resamples; ranks and residuals recomputed; percentile pointwise95CI; no selection adjustment','seed':20260920,'sensitivity':'exclude previously recorded warning IDs; effect-only comparison; no model selection','warning_ids':['FUSCCTNBC030','FUSCCTNBC044','FUSCCTNBC140'],'minimum_n':20,'stop_conditions':['missing variable','nonunique author patient ID','rank deficient design','constant residual'],'limitations':['same cohort is not independent validation','conditioning on other genes is not purity adjustment or causal mediation','data-selected followup P and q remain exploratory','not contribution fractions or proof of synergy'],'server_input':'/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A/20260920T090922Z_all174_external117_resources_v1/private/author_patient_joined_values.tsv'}
 (OUT/'analysis_spec.json').write_text(json.dumps(spec,indent=2)+'\n',encoding='utf-8')
 statcols=(ROOT/'templates/statistical_result.tsv').read_text().strip().split('\t')
 pending=[]
 for m in models:
  r=dict.fromkeys(statcols,'NA');r.update(cancer='BRCA',cohort='FUSCC_TNBC',stage_id='04_ROBUSTNESS',run_id=RUN,analysis_version='glutamine_context_v1',analysis_type='planned_partial_rank_'+m['model'],metabolite_key='KEGG:C00064',metabolite_name='glutamine',gene=m['gene'],unit='author_patient_id',effect_type='partial_rank_rho',test_family='exploratory_conditional6',status='ACCESS_BLOCKED',reason='server165_SSH_connection_timeout_three_attempts;not_computed',source_id='FUSCC_BRCA_2022',covariates=';'.join(m['covariates']))
  pending.append(r)
 write(OUT/'conditional_analysis_status.tsv',statcols+['covariates'],pending)
 manifest=[{'source':str(x.relative_to(ROOT)),'sha256':sha(x)} for x in [OLD/'external_all174_with_CAMP_context.tsv',ROOT/'code/brca_glutamine_followup_receipt_v1.py']]
 meta=OUT/'public_resource_metadata';meta.mkdir(exist_ok=True)
 url='https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA708216&result=read_run&fields=run_accession,sample_accession,experiment_accession,sample_title,library_strategy&format=tsv'
 p=meta/'PRJNA708216_runs.tsv'
 if not p.exists():
  x=requests.get(url,timeout=25);x.raise_for_status();p.write_bytes(x.content)
 manifest.append(dict(source=url,sha256=sha(p)))
 _,rr=read(p);assert len(rr)==2
 for a in ['SRX10312286','SRX10312287','SAMN18221691']:
  u='https://www.ebi.ac.uk/ena/browser/api/xml/'+a;f=meta/(a+'.xml')
  if not f.exists():
   x=requests.get(u,timeout=25);x.raise_for_status();f.write_bytes(x.content)
  manifest.append(dict(source=u,sha256=sha(f)))
 resource=[]
 for r in rr:
  t=ET.fromstring((meta/(r['experiment_accession']+'.xml')).read_bytes())
  title=t.findtext('.//EXPERIMENT/TITLE')
  resource.append(dict(gene='SLC6A8',accession='PRJNA708216',run=r['run_accession'],experiment=r['experiment_accession'],sample=r['sample_accession'],experiment_title=title,model='MDA-MB-231',perturbation='oxygen condition, not gene-specific',run_to_condition='UNRESOLVED',biological_replicates='UNVERIFIED',status='NEEDS_REVIEW',decision='Do not assign Sample_A/B to oxygen by order; do not infer independent repeats from runs or shared BioSample; no raw reprocessing'))
 write(OUT/'SLC6A8_resource_design.tsv',list(resource[0]),resource)
 write(OUT/'source_manifest.tsv',['source','sha256'],manifest)
 val={'status':'PARTIAL','reused_primary_significant_relations':len(sig),'opposite_point_estimate_relations':len(opposite),'all_opposite_CAMP_q_ge005':True,'same_sign_dual_original_FDR_genes':['ASNS','NNMT'],'new_patient_statistics':0,'new_conditional_rows':'ACCESS_BLOCKED','ENA_runs':len(rr),'ENA_unique_BioSamples':len({r['sample_accession'] for r in rr}),'condition_assignment_verified':False,'patient_data_downloaded':False,'server_archive_synced':False}
 (OUT/'validation.json').write_text(json.dumps(val,indent=2)+'\n',encoding='utf-8')
 sf,sr=read(ROOT/'coordination/stages/BRCA.tsv');sr=[r for r in sr if r['run_id']!=RUN]
 for stage,status,scope in [('04_ROBUSTNESS','ACCESS_BLOCKED','Conditional6 plan frozen; server SSH timeout; no new patient statistics'),('05_FUNCTION','PARTIAL','SLC6A8 ENA two runs one BioSample; condition labels and biological replication unresolved'),('07_INTEGRATION','DONE','Reuse original statistics; distinguish dual FDR from sign agreement and nonsignificant opposite estimates')]:
  vals=['BRCA',stage,RUN,'glutamine_context_v1',status,scope,OUT.relative_to(ROOT).as_posix(),'code/brca_glutamine_followup_receipt_v1.py','analysis/brca-functional-review-20260919','Post-selection exploratory; original q unchanged; server unavailable','Resume fixed conditional analysis when server reachable; resolve resource labels without guessing']
  sr.append(dict(zip(sf,vals)))
 write(ROOT/'coordination/stages/BRCA.tsv',sf,sorted(sr,key=lambda r:(r['stage_id'],r['run_id'])))
 (OUT/'.gitattributes').write_text('* -text\n',encoding='utf-8')
 for f in OUT.iterdir():
  if f.is_file():f.write_bytes(f.read_bytes().replace(b'\r\n',b'\n'))
 checksum=OUT/'checksums.sha256'
 checksum.write_bytes(''.join(f'{sha(p)}  {p.relative_to(OUT).as_posix()}\n' for p in sorted(OUT.rglob('*')) if p.is_file() and p!=checksum).encode('utf-8'))
 print(json.dumps(val))
if __name__=='__main__':main()
