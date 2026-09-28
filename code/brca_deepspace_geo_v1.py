"""Use original GEO archive for remaining single-specimen datasets."""
import re,json,hashlib
from brca_deepspace_fetch_v1 import D,get
from brca_gse210616_stream_to_server import remote,transfer
for sid,gsm in [('DSID000593','GSM7757981'),('DSID000595','GSM7757983')]:
 d=D+'/'+sid;remote("mkdir -p '"+d+"/spatial'")
 r=get(f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=quick');transfer([r.content],d+'/metadata.soft')
 records=[]
 for u in re.findall(r'!Sample_supplementary_file(?:_\d+)? = (\S+)',r.text):
  match=next((s for s in ['barcodes.tsv.gz','features.tsv.gz','matrix.mtx.gz','tissue_hires_image.png.gz','scalefactors_json.json.gz','tissue_positions_list.csv.gz'] if u.endswith(s)),None)
  if not match:continue
  u=u.replace('ftp://','https://');r=get(u);name=match if match in ['barcodes.tsv.gz','features.tsv.gz','matrix.mtx.gz'] else 'spatial/'+match
  transfer([r.content],d+'/'+name);records.append({'file':name,'url':u,'sha256':hashlib.sha256(r.content).hexdigest(),'bytes':len(r.content)})
 transfer([json.dumps({'id':sid,'files':records}).encode()],d+'/manifest.json');print(sid,'GEO_DONE',flush=True)
