"""Pinned Ensembl110 orthology; exact mouse symbol joins; preserve all117."""
from pathlib import Path
import argparse,json,requests,time,concurrent.futures,hashlib
import pandas as pd,numpy as np
BASE='https://e110.rest.ensembl.org'
def save(x,p):x.to_csv(p,sep='\t',index=False,na_rep='NA',lineterminator='\n')
def main(root):
 src=root/'source';pub=root/'public';cache=src/'ensembl110';cache.mkdir(exist_ok=True)
 genes=json.loads((src/'genes117.json').read_text());assert len(genes)==len(set(genes))==117
 def get(url,path):
  if path.exists():return json.loads(path.read_text())
  for attempt in range(3):
   try:
    x=requests.get(url,headers={'Accept':'application/json'},timeout=25)
    if x.status_code==429:time.sleep(min(float(x.headers.get('Retry-After',2)),20));continue
    x.raise_for_status();j=x.json();path.write_text(json.dumps(j));time.sleep(.12);return j
   except requests.RequestException:
    if attempt==2:raise
    time.sleep(1)
  raise RuntimeError('rate_limit_unresolved')
 release=get(BASE+'/info/data?content-type=application/json',cache/'release.json');assert release['releases']==[110]
 def one(g):
  row=dict(human_gene=g,human_ensembl='NA',mouse_ensembl='NA',mouse_gene='NA',homology_type='NA',status='NEEDS_REVIEW',reason='NA',source_release=110)
  try:
   j=get(BASE+'/homology/symbol/homo_sapiens/'+g+'?target_species=mus_musculus;type=orthologues;content-type=application/json',cache/(g+'_homology.json'));data=j.get('data',[])
   if len(data)!=1:
    exact=[]
    for d in data:
     lookup=get(BASE+'/lookup/id/'+d['id']+'?content-type=application/json',cache/(d['id']+'_human_lookup.json'))
     if lookup.get('species')=='homo_sapiens' and lookup.get('object_type')=='Gene' and lookup.get('display_name')==g:exact.append(d)
    data=exact
    if len(data)!=1:row['reason']='nonunique_or_absent_exact_human_display_name';return row
   row['human_ensembl']=data[0]['id'];h=data[0]['homologies'];row['n_mouse_homologies']=len(h)
   if len(h)!=1 or h[0]['type']!='ortholog_one2one':row['reason']='not_unique_one2one';row['mouse_ensembl']=';'.join(a['target']['id'] for a in h) or 'NA';return row
   t=h[0]['target'];assert t['species']=='mus_musculus';row['mouse_ensembl']=t['id'];row['homology_type']=h[0]['type']
   k=get(BASE+'/lookup/id/'+t['id']+'?content-type=application/json',cache/(t['id']+'_lookup.json'));assert k['species']=='mus_musculus' and k['object_type']=='Gene'
   row['mouse_gene']=k.get('display_name','NA');row['status']='DONE' if row['mouse_gene']!='NA' else 'NEEDS_REVIEW';row['reason']='Ensembl110_one2one;mouse_symbol_from_target_ID_lookup;not_case_conversion'
  except Exception as e:row.update(status='ACCESS_BLOCKED',reason=type(e).__name__+':'+str(e)[:150])
  return row
 rows=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
  for i,r in enumerate(ex.map(one,genes)):
   rows.append(r)
   if (i+1)%20==0:save(pd.DataFrame(rows),pub/'orthology_checkpoint.tsv');print('mapped',i+1,flush=True)
 ortho=pd.DataFrame(rows);save(ortho,pub/'orthology117.tsv')
 old=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A/20260921T025937Z_asns_gls_intervention_v1/public')
 coverage=pd.read_csv(old/'all_gene_test_coverage.tsv',sep='\t');oldtables={}
 for s in ['Tumor','Lung']:
  for c in ['shAsns1','shAsns2']:oldtables[(s,c)]=pd.read_csv(old/(s+'_'+c+'_all_DE.tsv'),sep='\t').set_index('gene')
 long=[];summ=[]
 for o in rows:
  for site in ['Tumor','Lung']:
   state=[]
   for con in ['shAsns1','shAsns2']:
    r=dict(human_gene=o['human_gene'],mouse_gene=o['mouse_gene'],site=site,construct=con,logFC=np.nan,p_value=np.nan,original_q_global4=np.nan,status=o['status'],reason=o['reason'],interpretation='response_to_Asns_intervention;not_self_gene_function_validation')
    if o['status']=='DONE':
     t=oldtables[(site,con)];g=o['mouse_gene']
     if g in t.index:
      x=t.loc[g];r.update(logFC=x.logFC,p_value=x.PValue,original_q_global4=x.q_global4,status='DONE',reason='exact_mouse_symbol_match;original49246_family_q_reused')
     else:r.update(status='NOT_EVALUABLE',reason='low_expression_not_tested' if ((coverage.gene==g)&(coverage.site==site)).any() else 'mouse_symbol_absent_in_author_matrix;no_alias_guess')
    long.append(r);state.append(r)
   if not all(x['status']=='DONE' for x in state):label='NOT_EVALUABLE_OR_MAPPING_UNRESOLVED'
   elif state[0]['logFC']*state[1]['logFC']<0:label='OPPOSITE_POINT_ESTIMATE_DIRECTIONS'
   elif all(x['original_q_global4']<.05 for x in state):label='BOTH_SUPPORTED_UP' if state[0]['logFC']>0 else 'BOTH_SUPPORTED_DOWN'
   elif any(x['original_q_global4']<.05 for x in state):label='SAME_DIRECTION_ONE_SUPPORTED'
   else:label='SAME_DIRECTION_NOT_JOINTLY_SUPPORTED'
   summ.append(dict(human_gene=o['human_gene'],mouse_gene=o['mouse_gene'],site=site,response_state=label,interpretation='Asns_intervention_response_only;not_equivalence_or_self_target_validation'))
 save(pd.DataFrame(long),pub/'all117_four_contrast_response.tsv');save(pd.DataFrame(summ),pub/'all117_two_construct_summary.tsv')
 manifest=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(cache.glob('*.json'))]+[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(old.glob('*_all_DE.tsv'))]
 save(pd.DataFrame(manifest),pub/'orthology_source_manifest.tsv')
 (pub/'mapping_validation.json').write_text(json.dumps(dict(genes=117,long_rows=len(long),site_summary_rows=len(summ),mapping_status=ortho.status.value_counts().to_dict(),no_new_DE_tests=True,old_FDR_unchanged=True),indent=2)+'\n')
 print((pub/'mapping_validation.json').read_text(),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args().root)
