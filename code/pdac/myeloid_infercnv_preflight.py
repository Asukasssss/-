"""Server-only input readiness. Raw metadata and per-unit tables remain private."""
from pathlib import Path
import csv,gzip,hashlib,json,re,sys
import pandas as pd
RUN=Path(__file__).resolve().parent
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
DATA=ROOT/'data/candidates'
assert RUN.parent==ROOT/'results/collaborative/PDAC/B'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
inputs=[DATA/'gse263733_pdac_scRNA_v0.1/GSE263733_Cell_annotation.txt.gz',DATA/'gse278688_pdac_scRNA_v0.1/GSE278688_sc_metadata.csv.gz',DATA/'GSE242230/GSE242230_annotations.txt.gz']
summ=[]
for c,p in zip(['GSE263733','GSE278688','GSE242230'],inputs):
    d=pd.read_csv(p,sep=',' if c=='GSE278688' else '\t');raw=d.copy()
    if c=='GSE263733':
        d=d[d.Origins.isin(['Pm0','Pm1'])].rename(columns={'Barcodes':'cell','Patients':'unit','Celltypes':'label'})
        d['role']=d.label.map({'Myeloid':'Myeloid','T cell':'T_ref','B cell':'B_ref'});d['subtype']=d.label
    elif c=='GSE278688':
        d=d[d.tissue=='Tumor'].rename(columns={d.columns[0]:'cell','patients':'unit','all_celltype':'label'})
        d['role']=d.label.map({'Myeloid':'Myeloid','CD4T':'T_ref','CD8T':'T_ref','B':'B_ref'});d['subtype']=d.label
    else:
        d=d[(d.dataset=='eusfnb')&(d.tissue_type=='tumor_primary')&(~d.filtered)&d.cell_type.notna()].rename(columns={'cell_id':'cell','sample_id':'unit','cell_type':'label','cell_type_specific':'subtype'})
        d['role']=d.label.map({'Monocyte/DC':'Myeloid','B cell':'B_ref'})
        d.loc[d.subtype.str.strip().isin(['CD4 T cell','CD8 T cell']),'role']='T_ref'
        d.loc[d.subtype.isin(['Malignant - Basal','Malignant - Classical']),'role']='Author_malignant_control'
    if c=='GSE263733':
        n=raw[(raw.Origins=='Pn')&(raw.Celltypes=='Myeloid')].rename(columns={'Barcodes':'cell','Patients':'unit','Celltypes':'label'})
        n['role']='Myeloid_normal_control';n['subtype']=n.label;d=pd.concat([d,n],ignore_index=True)
    elif c=='GSE278688':
        n=raw[(raw.tissue=='Adjacent_normal')&(raw.all_celltype=='Myeloid')].rename(columns={raw.columns[0]:'cell','patients':'unit','all_celltype':'label'})
        n['role']='Myeloid_normal_control';n['subtype']=n.label;d=pd.concat([d,n],ignore_index=True)
    assert d.cell.is_unique and d[['cell','unit']].notna().all().all()
    d=d[d.role.notna()].copy()
    d.to_csv(RUN/'private'/f'{c}_metadata.tsv',sep='\t',index=False)
    counts=pd.crosstab(d.unit,d.role).reindex(columns=['Myeloid','T_ref','B_ref','Author_malignant_control'],fill_value=0)
    counts['eligible']=(counts.Myeloid>=50)&(counts.T_ref>=50)&(counts.B_ref>=50)
    counts.to_csv(RUN/'private'/f'{c}_unit_counts.tsv',sep='\t')
    summ.append({'cancer':'PDAC','cohort':c,'stage_id':'06_EXTERNAL','run_id':RUN.name,'analysis_version':'myeloid_infercnv_v1','n_units_total':len(counts),'n_units_reference_eligible':int(counts.eligible.sum()),'n_myeloid_total':int(counts.Myeloid.sum()),'n_myeloid_reference_eligible':int(counts.loc[counts.eligible,'Myeloid'].sum()),'status':'DONE','reason':'Metadata preflight only;CNV not yet run;at least50 myeloid,T and B cells per author unit'})
pd.DataFrame(summ).to_csv(RUN/'public/reference_readiness.tsv',sep='\t',index=False)
gtf=Path('/public3/xuzx/Cancer/lung_cancer_pack/raw/bulk_transcriptomics/xena_gdc_tcga/gencode.v36.annotation.gtf.gz')
genes=[]
with gzip.open(gtf,'rt') as f:
    for line in f:
        if line.startswith('#'):continue
        a=line.rstrip().split('\t')
        if a[2]!='gene' or a[0] not in {'chr'+str(i) for i in range(1,23)}:continue
        name=re.search(r'gene_name "([^"]+)"',a[8])
        if name:genes.append((name[1],a[0],int(a[3]),int(a[4])))
g=pd.DataFrame(genes,columns=['gene','chr','start','end'])
slc=g[g.gene=='SLC6A6'];assert len(slc)==1 and slc.iloc[0]['chr']=='chr3'
g=g[~g.gene.duplicated(keep=False)&(g.chr!='chr3')].copy()
g['order']=g.chr.str[3:].astype(int);g=g.sort_values(['order','start']).drop(columns='order')
g.to_csv(RUN/'private/gene_order_no_chr3.tsv',sep='\t',index=False,header=False)
pd.DataFrame([{'path':str(p),'sha256':sha(p)} for p in inputs+[gtf]]).to_csv(RUN/'public/source_manifest.tsv',sep='\t',index=False)
(RUN/'public/preflight_validation.json').write_text(json.dumps({'status':'PASS','unique_cell_keys':True,'gene_order_assembly':'GRCh38 GENCODE v36','gene_order_count':len(g),'SLC6A6_chr':'chr3','excluded_chr3':True,'new_CNV_results':False,'clinical_donor_identity_recertified':False},indent=2))
print(json.dumps(summ),flush=True)
