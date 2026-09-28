"""Finalize source provenance without repeating numerical analyses."""
import sys,pathlib,json,hashlib
import pandas as pd
r=pathlib.Path(sys.argv[1]);old=pathlib.Path(sys.argv[2]);commit=sys.argv[3]
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
base='https://cells.ucsc.edu/multiomic-gbm/'
urls={'SNUH_collection.json':base+'dataset.json','SNUH_scrna_dataset.json':base+'scrna/dataset.json',
 'SNUH_exprMatrix.json':base+'scrna/exprMatrix.json','SNUH_meta_proxy.tsv':base+'scrna/meta.tsv',
 'SNUH_barcodes.tsv.gz':base+'scrna/barcodes.tsv.gz','LYPLA1.range.zlib':base+'scrna/exprMatrix.bin;bytes=744376841-744581768',
 'SNUH_clinical_supplement.xlsx':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-69716-2/MediaObjects/41467_2026_69716_MOESM5_ESM.xlsx',
 'SS2_metadata.xlsx':'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8452nnn/GSM8452708/suppl/GSM8452708_All_cells_metadata.xlsx',
 'SS2_expression.rds.gz':'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8452nnn/GSM8452708/suppl/GSM8452708_PvR_GBM_SS2_SuvaLab.rds.gz'}
paths=[p for p in (r/'source').iterdir() if p.is_file() and not p.name.startswith('meta.part.') and p.name not in ['SNUH_meta.tsv','SNUH_meta_parallel.tsv']]
paths+=list((r/'code').glob('gbm_lypla1*'))+[old/'private/sc_donor_profiles_private.tsv']
df=pd.DataFrame([dict(source_path=str(p),source_url=urls.get(p.name,'derived_local_or_repository_code;see_analysis_spec'),bytes=p.stat().st_size,sha256=sha(p),scope='server_only_input_or_code') for p in paths])
df.to_csv(r/'public/source_manifest.tsv',sep='\t',index=False)
spec=json.loads((r/'public/analysis_spec.json').read_text());spec['input_code_base_commit']=commit
spec['source_reference']='SNUH paper DOI10.1038/s41467-026-69716-2;CARE DOI10.1038/s41588-025-02168-4'
spec['source_snapshot']='portal byte hashes pinned in source_manifest.tsv;author CARE annotation commit383b0a320f0993a4529167cdb591dac63e103674'
(r/'public/analysis_spec.json').write_text(json.dumps(spec,indent=2))
print('Source and code hashes finalized',len(df))
