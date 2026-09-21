"""Resumable official GEO processed sample records and TCGA RNA via cBioPortal."""
from pathlib import Path
import requests,sys,json,time,concurrent.futures,re
import pandas as pd
r=Path(sys.argv[1]);src=r/'source'
def get(url,path,params=None):
 if path.exists():return path.read_text()
 for i in range(3):
  try:
   x=requests.get(url,params=params,timeout=50);x.raise_for_status();path.write_text(x.text);return x.text
  except Exception:
   if i==2:raise
   time.sleep(2)
def geo(acc,view='full',targ='self'):
 return get('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi',src/(acc+'_'+targ+'_'+view+'.soft'),dict(acc=acc,targ=targ,form='text',view=view))
if sys.argv[2]=='geo':
 geo('GPL570');geo('GSE58212','brief','all')
 text=(src/'GSE42568_series.soft').read_text();ids=[l.split(' = ')[-1] for l in text.splitlines() if l.startswith('!Series_sample_id')]
 def one(a):
  try:geo(a);return a,'DONE'
  except Exception as e:return a,str(e)
 rows=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
  for i,x in enumerate(ex.map(one,ids)):
   rows.append(x)
   if (i+1)%20==0:print('GSE42568 processed samples',i+1,flush=True)
 (r/'public/geo_download_status.json').write_text(json.dumps(dict(samples_requested=len(rows),samples_downloaded=sum(x[1]=='DONE' for x in rows),failures=sum(x[1]!='DONE' for x in rows)),indent=2))
else:
 api='https://www.cbioportal.org/api';study='brca_tcga_pub2015';profile=study+'_rna_seq_v2_mrna'
 get(api+'/studies/'+study+'/molecular-profiles',src/'tcga_profiles.json')
 samples=json.loads(get(api+'/studies/'+study+'/samples',src/'tcga_samples.json',dict(projection='DETAILED',pageSize=10000)))
 d=pd.read_excel(src/'Tang_S2.xlsx',header=None);labels=[str(x).strip() for x in d.iloc[0,3:]];patients=['TCGA-'+x for x in labels if re.fullmatch(r'[A-Z0-9]{2}-[A-Z0-9]{4}',x)];assert len(patients)==len(set(patients))==25
 selected=[x for x in samples if x['patientId'] in patients and x['sampleType']=='Primary Solid Tumor'];assert len({x['patientId'] for x in selected})==len(selected)
 (src/'Tang_tcga_sample_join.json').write_text(json.dumps(selected))
 rel=pd.read_csv(src/'relations174.tsv',sep='\t');genes=sorted(rel.gene.unique());info=[]
 def gene(g):
  try:
   j=json.loads(get(api+'/genes/'+g,src/('tcga_gene_'+g+'.json')))
   if j.get('hugoGeneSymbol')==g:return j
  except Exception:pass
  return None
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
  for j in ex.map(gene,genes):
   if j:info.append(j)
 for i in range(0,len(info),20):
  gs=info[i:i+20];p=src/('tcga_rna_'+str(i)+'.json')
  if not p.exists():
   x=requests.post(api+'/molecular-profiles/'+profile+'/molecular-data/fetch',params={'projection':'DETAILED'},json={'entrezGeneIds':[x['entrezGeneId'] for x in gs],'sampleIds':[x['sampleId'] for x in selected]},timeout=60);x.raise_for_status();p.write_text(x.text)
  print('TCGA RNA batch',i,flush=True)
 (r/'public/tcga_download_summary.json').write_text(json.dumps(dict(patient_barcodes_in_metabolomics=25,unique_primary_samples_matched=len(selected),gene_symbols_resolved=len(info),profile=profile),indent=2))
