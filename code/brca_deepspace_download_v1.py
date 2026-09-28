"""Download expression only after expert-label eligibility, without local data files."""
import json,hashlib,time,concurrent.futures,bs4
from brca_deepspace_fetch_v1 import B,D,get
from brca_gse210616_stream_to_server import transfer,remote
def one(rec):
 x=rec['metadata'];sid=x['ID'];d=D+'/'+sid;remote("mkdir -p '"+d+"/spatial'")
 page=bs4.BeautifulSoup(get(B+'/database/'+sid).text,'html.parser')
 available=[a['href'].split('/'+sid+'/')[-1] for a in page.find_all('a',href=True) if '/download/'+sid+'/' in a['href']]
 files=[f for f in available if f in ['filtered_feature_bc_matrix.h5','barcodes.tsv.gz','features.tsv.gz','matrix.mtx.gz','spatial/tissue_hires_image.png','spatial/scalefactors_json.json','spatial/tissue_positions_list.csv']]
 out=[]
 for f in files:
  u=B+'/download/'+sid+'/'+f;r=get(u);transfer([r.content],d+'/'+f);out.append({'file':f,'url':u,'sha256':hashlib.sha256(r.content).hexdigest(),'bytes':len(r.content)})
 gsm=x['sample_id'];u=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=quick';r=get(u);transfer([r.content],d+'/metadata.soft')
 u=B+'/api/sample/'+sid+'/spatial-data';r=get(u);transfer([r.content],d+'/spatial_api.json')
 result={'id':sid,'files':out};transfer([json.dumps(result).encode()],d+'/manifest.json')
 print(sid,gsm,'DONE',flush=True);return result
if __name__=='__main__':
 j=json.loads(remote("cat '"+D+"/expert_annotations.json'"))
 rec=[x for x in j['records'] if 'Normal duct' in x['labels'] and x['metadata']['series_id'] in ['GSE210616','GSE242311']]
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:r=list(ex.map(one,rec))
 transfer([json.dumps(r,indent=2).encode()],D+'/download_manifest.json')
