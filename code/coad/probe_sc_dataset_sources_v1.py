"""Read-only source availability probe. No expression/cell/patient payload downloads."""
import json,urllib.request
j=json.load(urllib.request.urlopen('https://zenodo.org/api/records/16631519',timeout=30))
for x in j['files']:
 print(json.dumps({'name':x['key'],'bytes':x['size'],'checksum':x['checksum'],'url':x['links']['self']}))
u='https://zenodo.org/api/records/16631519/files/MUI_Innsbruck-adata.h5ad/content'
with urllib.request.urlopen(urllib.request.Request(u,method='HEAD'),timeout=30) as f:print('HEAD',f.status,f.headers.get('Content-Length'))
