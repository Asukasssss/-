"""Retrieve only processed counts and explicit author metadata on server165."""
from pathlib import Path
import argparse,re,json,requests,hashlib,time,pandas as pd

def main(root):
 src=root/'source';pub=root/'public';src.mkdir(exist_ok=True);pub.mkdir(exist_ok=True)
 queries=['(ASNS[All Fields] OR "asparagine synthetase"[All Fields]) AND gse[Entry Type]','(GLS[All Fields] OR GLS1[All Fields] OR glutaminase[All Fields]) AND (breast[All Fields] OR mammary[All Fields]) AND gse[Entry Type]']
 for i,query in enumerate(queries):
  sp=src/f'search{i}.json';su=src/f'summary{i}.json'
  if not sp.exists():
   x=requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',params={'db':'gds','term':query,'retmax':100,'retmode':'json','tool':'CAMP_BRCA_resource_lookup'},timeout=60);x.raise_for_status();sp.write_text(json.dumps(x.json()));time.sleep(.4)
  j=json.loads(sp.read_text())['esearchresult'];assert int(j['count'])<=100, 'Increase recorded retrieval bound before claiming complete search results'
  if not su.exists():
   x=requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi',params={'db':'gds','id':','.join(j['idlist']),'retmode':'json'},timeout=60);x.raise_for_status();su.write_text(json.dumps(x.json()));time.sleep(.4)
 manifest=[];qcs=[]
 for acc in ['GSE104967','GSE104966','GSE104964']:
  p=src/(acc+'.txt');url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+acc+'&targ=all&form=text&view=full'
  if not p.exists():
   x=requests.get(url,timeout=60);x.raise_for_status();p.write_bytes(x.content)
  text=p.read_text();assert '^SERIES = '+acc in text
  u=re.search(r'!Series_supplementary_file = (\S+)',text).group(1).replace('ftp://','https://');cp=src/(acc+'_counts.txt.gz')
  if not cp.exists():
   x=requests.get(u,timeout=120);x.raise_for_status();cp.write_bytes(x.content)
  a=pd.read_csv(cp,sep='\t',index_col=0);assert a.index.is_unique
  for col in a.columns:qcs.append(dict(dataset=acc,column=col,library_total=int(a[col].sum()),nonzero_genes=int((a[col]>0).sum()),status='NEEDS_REVIEW' if a[col].sum()==0 else 'DONE',reason='all_zero_library' if a[col].sum()==0 else 'nonempty_count_library_not_full_QC'))
  for f,url2 in [(p,url),(cp,u)]:manifest.append(dict(path=str(f),url=url2,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
 text=(src/'GSE104966.txt').read_text();assert 'Counts per gene for each biological replicate' in text
 desc=re.search(r'Supplementary_files_format_and_content: Counts per gene for each biological replicate\. Columns correspond to the following samples: (.*)',text).group(1)
 entries=re.split(r'\d+\)',desc)[1:];assert len(entries)==24
 order=[]
 samples={}
 for block in text.split('^SAMPLE = ')[1:]:
  gsm=block.splitlines()[0].strip();title=re.search(r'!Sample_title = (.*)',block).group(1).strip();key=re.sub(r'[^a-z0-9]','',title.lower().replace('replicate',''));samples[key]=gsm
 for i,e in enumerate(entries):
  e=e.strip();site='Tumor' if i<12 else 'Lung';group=['control','shAsns1','shAsns2'][(i%12)//4];rep='ABCD'[i%4]
  assert e.startswith(site) and ('Renilla' in e if group=='control' else 'ASNS' in e)
  # Published 21st description lacks A; exact GEO sample title explicitly supplies it.
  title=site+' '+('shRenilla' if group=='control' else 'shAsns-'+group[-1])+' Replicate '+rep
  key=re.sub(r'[^a-z0-9]','',title.lower().replace('replicate',''));assert key in samples
  stated=re.sub(r'[^a-z0-9]','',e.lower())
  assert stated==key or (i==20 and stated==key[:-1]), (i,e,title)
  order.append(dict(column=i+1,site=site,group=group,replicate_label=rep,author_sample=samples[key],author_column_description=e,assignment_basis='explicit_author_numbered_column_description_and_GEO_title;not_inferred_pairing'))
 pd.DataFrame(order).to_csv(src/'author_column_assignment.tsv',sep='\t',index=False)
 pd.DataFrame(qcs).to_csv(pub/'processed_count_coverage.tsv',sep='\t',index=False)
 pd.DataFrame(manifest).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
 spec=dict(version='asns_intervention_v1',primary_dataset='GSE104966',species='Mus musculus',cell_model='4T1-T;recovered_from_tumor_or_lung;6TG_selected_recultured',unit='author_labeled_biological_replicate',n_per_group=4,sites=['Tumor','Lung'],constructs=['shAsns1','shAsns2'],control='shRenilla',pairing='none;shared letters do not imply paired samples',analysis='edgeR QL F;each source separate model~0+group;robust dispersion and QL;TMM;filterByExpr min.count10 min.total.count15',family='all_tested_genes_x_four_contrasts_joint_BH;within_contrast_q_also_retained',no_outcome_based_exclusion=True,focus_genes=['Asns','Gls','Gls2','Glul'],frozen_before_DE=True,excluded_dataset={'GSE104967':'all-zero control Sample2;only one effective control;no DE'},other_dataset={'GSE104964':'Asns-silenced cells plus/minus asparagine;nutrient contrast not new gene knockdown;not_run_this_batch'},limitations=['author replicate independence not independently verified by mouse IDs','reculture and 6TG selection limit tissue inference','RNA only;not metabolite flux or human TNBC validation','two constructs share controls and study;not two independent studies'],source_commit='d328581')
 sp=root/'analysis_spec.json'
 if sp.exists():assert json.loads(sp.read_text())==spec
 else:sp.write_text(json.dumps(spec,indent=2)+'\n')
 print('Prepared explicit author columns and frozen specification')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args().root)
