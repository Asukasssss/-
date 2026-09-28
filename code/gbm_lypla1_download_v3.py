"""Optional reproducible acquisition. Run on server165 only; no target selection."""
import pathlib,sys,requests,shutil,gzip
ROOT=pathlib.Path(sys.argv[1]);s=ROOT/'source';s.mkdir(exist_ok=True)
BASE='https://cells.ucsc.edu/multiomic-gbm/'
URLS={
 'SNUH_collection.json':BASE+'dataset.json',
 'SNUH_scrna_dataset.json':BASE+'scrna/dataset.json',
 'SNUH_exprMatrix.json':BASE+'scrna/exprMatrix.json',
 'SNUH_meta_proxy.tsv':BASE+'scrna/meta.tsv',
 'SNUH_barcodes.tsv.gz':BASE+'scrna/barcodes.tsv.gz',
 'SNUH_clinical_supplement.xlsx':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-69716-2/MediaObjects/41467_2026_69716_MOESM5_ESM.xlsx',
 'SS2_metadata.xlsx':'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8452nnn/GSM8452708/suppl/GSM8452708_All_cells_metadata.xlsx',
 'SS2_expression.rds.gz':'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8452nnn/GSM8452708/suppl/GSM8452708_PvR_GBM_SS2_SuvaLab.rds.gz',
}
for name,url in URLS.items():
    if (s/name).exists():continue
    r=requests.get(url,stream=True,timeout=(20,90));r.raise_for_status()
    with (s/(name+'.partial')).open('wb') as f:
        for chunk in r.iter_content(1024*1024):f.write(chunk)
    (s/(name+'.partial')).rename(s/name)
if not (s/'SS2_expression.rds').exists():
    with gzip.open(s/'SS2_expression.rds.gz','rb') as i,(s/'SS2_expression.rds').open('wb') as o:shutil.copyfileobj(i,o)
# For original input identity, compare hashes with the published source_manifest.tsv.
# Rerun validation requires a second independent download of meta.tsv saved as SNUH_meta_parallel.tsv.
