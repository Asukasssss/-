"""Public HD 3-prime / Kinnex source data streamed only to server165."""
import requests,json,hashlib,datetime,sys
from brca_gse210616_stream_to_server import remote,transfer
B='https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/DATA-RevioSPRQ-Kinnex-VisiumHD-humanBreast/'
D='/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/BRCA_HD3_Kinnex'
RUN='20260928T040500Z_hd3_lypla1_v1'
OUT='/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/BRCA/A/'+RUN
if __name__=='__main__':
 remote(f"mkdir -p '{D}' '{OUT}'; mkdir '{OUT}/.running'")
 records=[]
 for name,u in [('genes.tsv',B+'5-CountMatrix/genes_seurat/genes.tsv'),('barcodes.tsv',B+'5-CountMatrix/genes_seurat/barcodes.tsv'),('matrix.mtx',B+'5-CountMatrix/genes_seurat/matrix.mtx'),('web_summary.html',B+'hBreast_web_summary.html'),('source_README.txt','https://downloads.pacbcloud.com/public/dataset/Kinnex-single-cell-RNA/README.txt')]:
  with requests.get(u,stream=True,timeout=(30,90)) as r:
   r.raise_for_status();h=hashlib.sha256();size=[0]
   def chunks():
    for b in r.iter_content(1024*1024):h.update(b);size[0]+=len(b);yield b
   transfer(chunks(),D+'/'+name)
   records.append({'file':name,'url':u,'sha256':h.hexdigest(),'bytes':size[0],'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
   print(name,size[0],'DONE',flush=True)
 transfer([json.dumps(records,indent=2).encode()],D+'/manifest.json')
