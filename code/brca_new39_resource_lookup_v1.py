"""Public full-text availability lookup; accession discovery is not dataset eligibility."""
from pathlib import Path
import sys,json,re,concurrent.futures
import xml.etree.ElementTree as ET
import pandas as pd,requests
R=Path(sys.argv[1]);(R/'fulltext').mkdir(exist_ok=True)
d=pd.concat([pd.read_csv(R/f,sep='\t',dtype=str) for f in ['literature_hits.tsv','literature_hits_focused.tsv']]).drop_duplicates(['gene','pmid','doi'])
pmids=['32366405','31581558','38953696','39700137','38078907','41271644','33224971','38646426','37296899','36909375','41992085','28065656','35351947','39515327','37253003','30205167','37974462','31164982','40370195','40821099']
d=d[d.pmid.isin(pmids)].drop_duplicates('pmid')
def one(r):
 o=dict(pmid=r.pmid,pmcid=r.pmcid,gene=r.gene,title=r.title,source_url=r.url,status='NOT_EVALUABLE',reason='no_PMC_identifier',accessions='',availability_text='',model_text='')
 if pd.isna(r.pmcid):return o
 p=R/'fulltext'/(r.pmcid+'.xml')
 try:
  if not p.exists():
   x=requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/'+r.pmcid+'/fullTextXML',timeout=35);x.raise_for_status();p.write_bytes(x.content)
  root=ET.fromstring(p.read_bytes());pars=[' '.join(e.itertext()) for e in root.iter('p')];o.update(status='DONE',reason='fulltext_retrieved;resource_eligibility_not_yet_verified',accessions=';'.join(sorted(set(re.findall(r'\b(?:GSE|PRJNA|E-MTAB-)\d+\b',' '.join(pars))))),availability_text=' | '.join(t for t in pars if re.search('GSE\d|data availability|data generated|data supporting|data that support|data are available|data is available|deposited',t,re.I))[:12000],model_text=' | '.join(t for t in pars if re.search('cell lines? (?:were|was|used)|MDA-MB|MCF-?7|4T1|NCI-H460|MDA231|MDA468',t,re.I))[:14000])
 except Exception as e:o.update(status='ACCESS_BLOCKED',reason=type(e).__name__+':'+str(e)[:120])
 return o
out=list(concurrent.futures.ThreadPoolExecutor(max_workers=3).map(one,d.itertuples()));pd.DataFrame(out).to_csv(R/'public_resource_lookup.tsv',sep='\t',index=False);print(pd.DataFrame(out)[['gene','pmid','status','accessions']].to_string(index=False))
# Official series-level metadata only, not the count matrix or any patient data.
gp=R/'sources/GSE283282_metadata.txt'
if not gp.exists():
 x=requests.get('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi',params={'acc':'GSE283282','targ':'self','form':'text','view':'full'},timeout=35);x.raise_for_status();assert '!Series_geo_accession = GSE283282' in x.text;gp.write_text(x.text,encoding='utf-8')
