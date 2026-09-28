"""Server-only, patient-balanced descriptive gene sources in three original cohorts."""
import argparse,csv,gzip,hashlib,json,platform,traceback
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
DATA=ROOT/'data/candidates'
SOURCE_FILES=set()
def read_table(p,**kw):SOURCE_FILES.add(p);return pd.read_csv(p,**kw)
def lines(p):
 SOURCE_FILES.add(p)
 with gzip.open(p,'rt') as f:return [x.rstrip('\r\n') for x in f]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(p,data):p.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
def broad(label):
 return {'Ductal cell':'Epithelial/Ductal','Ductal':'Epithelial/Ductal','Malignant':'Epithelial/Ductal','Fibroblast':'Fibroblast/CAF','CAF':'Fibroblast/CAF','Stellate':'Fibroblast/CAF','T cell':'T/NK','CD4T':'T/NK','CD8T':'T/NK','NK':'T/NK','NK/T cell':'T/NK','B cell':'B/Plasma','B':'B/Plasma','Plasma':'B/Plasma','Monocyte/DC':'Myeloid','Endothelial cell':'Endothelial','Acinar cell':'Acinar','Platelets':'Platelet'}.get(label,label)
def sparse_extract(matrix,features,barcodes,meta,genes):
 SOURCE_FILES.add(matrix)
 assert len(set(barcodes))==len(barcodes) and meta.index.is_unique
 assert meta.index.isin(barcodes).all()
 pos=pd.Index(barcodes).get_indexer(meta.index)
 cmap=np.full(len(barcodes),-1,dtype=np.int32);cmap[pos]=np.arange(len(pos))
 multiplicity=Counter(features)
 gmap=np.full(len(features),-1,dtype=np.int32);present=[]
 lookup={g:i for i,g in enumerate(genes)}
 for i,g in enumerate(features):
  if g in lookup and multiplicity[g]==1:gmap[i]=lookup[g];present.append(g)
 target=np.zeros((len(genes),len(pos)),dtype=np.int64);lib=np.zeros(len(pos),dtype=np.int64)
 with gzip.open(matrix,'rt') as f:
  assert f.readline().startswith('%%MatrixMarket matrix coordinate integer')
  line=f.readline()
  while line.startswith('%'):line=f.readline()
  nr,nc,nnz=map(int,line.split());assert (nr,nc)==(len(features),len(barcodes));seen=0
  for block in pd.read_csv(f,sep=' ',header=None,names=['row','col','value'],dtype=np.int64,chunksize=3000000):
   a=block.to_numpy();seen+=len(a)
   assert (a[:,2]>=0).all() and a[:,0].min()>=1 and a[:,0].max()<=nr and a[:,1].min()>=1 and a[:,1].max()<=nc
   ci=cmap[a[:,1]-1];mask=ci>=0;a=a[mask];ci=ci[mask]
   np.add.at(lib,ci,a[:,2]);gi=gmap[a[:,0]-1];m=gi>=0
   np.add.at(target,(gi[m],ci[m]),a[m,2])
  assert seen==nnz
 assert np.all(target.sum(axis=0)<=lib) and (lib>0).all()
 return target,lib,set(present)
def dense263(genes):
 folder=DATA/'gse263733_pdac_scRNA_v0.1';ann=read_table(folder/'GSE263733_Cell_annotation.txt.gz',sep='\t')
 ann=ann[ann.Origins.isin(['Pm0','Pm1'])].copy().set_index('Barcodes');ann=ann.rename(columns={'Patients':'unit','Celltypes':'celltype'})
 assert ann.index.is_unique and ann[['unit','celltype']].notna().all().all()
 p=folder/'GSE263733_Raw_counts.txt.gz';SOURCE_FILES.add(p);target=np.zeros((len(genes),len(ann)),dtype=np.int64);lib=np.zeros(len(ann),dtype=np.int64);seen=set();present=set();lookup={g:i for i,g in enumerate(genes)}
 with gzip.open(p,'rt') as f:
  bars=f.readline().rstrip().split('\t');assert len(set(bars))==len(bars);idx=pd.Index(bars).get_indexer(ann.index);assert (idx>=0).all()
  for line in f:
   gene,values=line.split('\t',1);assert gene not in seen;seen.add(gene)
   x=np.fromstring(values,sep='\t',dtype=np.int64);assert len(x)==len(bars) and (x>=0).all();x=x[idx];lib+=x
   if gene in lookup:target[lookup[gene]]=x;present.add(gene)
 assert (lib>0).all() and np.all(target.sum(axis=0)<=lib)
 return ann,target,lib,present
def sparse278(genes):
 folder=DATA/'gse278688_pdac_scRNA_v0.1';ann=read_table(folder/'GSE278688_sc_metadata.csv.gz',index_col=0);ann=ann[ann.tissue=='Tumor'].copy().rename(columns={'patients':'unit','all_celltype':'celltype'})
 bars=lines(folder/'GSE278688_sc_barcodes.tsv.gz');features=[x.split('\t')[0] for x in lines(folder/'GSE278688_sc_features.tsv.gz')]
 a,lib,present=sparse_extract(folder/'GSE278688_sc_matrix.mtx.gz',features,bars,ann,genes)
 return ann,a,lib,present
def sparse242(genes):
 folder=DATA/'GSE242230';ann=read_table(folder/'GSE242230_annotations.txt.gz',sep='\t')
 ann=ann[(ann.dataset=='eusfnb')&(ann.tissue_type=='tumor_primary')&(~ann.filtered)&ann.cell_type.notna()].copy().set_index('cell_id').rename(columns={'sample_id':'unit','cell_type':'celltype'})
 aa=[];xx=[];ll=[];present_all=None
 for unit,meta in ann.groupby('unit',sort=True):
  matches=list((folder/'raw_matrices').glob('*_'+unit+'_matrix.mtx.gz'));assert len(matches)==1;matrix=matches[0];prefix=str(matrix)[:-len('matrix.mtx.gz')]
  bars=[unit+'_'+s for s in lines(Path(prefix+'barcodes.tsv.gz'))];features=[s.split('\t')[1] for s in lines(Path(prefix+'features.tsv.gz'))]
  # Ambiguous duplicate symbols are excluded from the target panel; all rows remain in total UMI.
  a,lib,present=sparse_extract(matrix,features,bars,meta,genes);aa.append(meta);xx.append(a);ll.append(lib)
  if present_all is None:present_all=present
  else:assert present_all==present
  print('GSE242230',unit,'done',flush=True)
 return pd.concat(aa),np.concatenate(xx,axis=1),np.concatenate(ll),present_all
def summarize(cohort,meta,target,lib,present,genes,run):
 assert meta.index.is_unique and meta[['unit','celltype']].notna().all().all()
 meta=meta.copy();meta['broad']=meta.celltype.map(broad);meta['position']=np.arange(len(meta));private=[];output=[]
 for scheme,column in [('author','celltype'),('broad','broad')]:
  rows=[]
  for (unit,ct),m in meta.groupby(['unit',column],sort=True):
   idx=m.position.to_numpy();counts=target[:,idx];nc=len(idx);total=int(lib[idx].sum());summed=counts.sum(axis=1);fraction=(counts>0).mean(axis=1)
   for i,gene in enumerate(genes):
    rows.append({'cohort':cohort,'scheme':scheme,'unit':unit,'celltype':ct,'gene':gene,'n_cells':nc,'library_umi':total,'sum_counts':int(summed[i]) if gene in present else np.nan,'log1p_cpm':float(np.log1p(summed[i]*1e6/total)) if gene in present else np.nan,'positive_fraction':float(fraction[i]) if gene in present else np.nan})
  pb=pd.DataFrame(rows);private.append(pb)
  for (gene,ct),g in pb.groupby(['gene','celltype'],sort=True):
   eligible=g[g.n_cells>=20];n=len(eligible);det=int((eligible.sum_counts>0).sum())
   status='GENE_ABSENT_OR_AMBIGUOUS' if gene not in present else 'INSUFFICIENT_UNITS' if n<3 else 'LOW_DETECTION' if det<3 else 'EVALUABLE'
   output.append({'cohort':cohort,'annotation_level':scheme,'gene':gene,'celltype':ct,'n_units_total':len(g),'n_units_ge20cells':n,'n_units_detected':det,'n_cells_total':int(g.n_cells.sum()),'mean_unit_log1p_cpm':float(eligible.log1p_cpm.mean()) if n and gene in present else np.nan,'median_unit_log1p_cpm':float(eligible.log1p_cpm.median()) if n and gene in present else np.nan,'mean_unit_positive_fraction':float(eligible.positive_fraction.mean()) if n and gene in present else np.nan,'status':status,'source_rank':np.nan})
 pd.concat(private).to_csv(run/(cohort+'_private_unit_pseudobulk.tsv'),sep='\t',index=False,na_rep='NA')
 out=pd.DataFrame(output)
 for (scheme,gene),g in out.groupby(['annotation_level','gene']):
  idx=g.index[g.status=='EVALUABLE'];out.loc[idx,'source_rank']=out.loc[idx,'mean_unit_log1p_cpm'].rank(ascending=False,method='min')
 cells=meta.groupby('celltype').agg(n_cells=('unit','size'),n_units=('unit','nunique')).reset_index();cells.insert(0,'cohort',cohort)
 return out,cells,{'cohort':cohort,'cells_included':len(meta),'units_included':int(meta.unit.nunique()),'target_genes':len(genes),'genes_in_matrix':len(present),'author_celltypes':int(meta.celltype.nunique()),'barcode_join':'PASS','raw_count_integrity':'PASS','missing_genes':sorted(set(genes)-present)}
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--code-commit',required=True);args=parser.parse_args()
 run=Path(__file__).resolve().parent;assert run.parent==ROOT/'results/collaborative/PDAC/B'
 with (run/'.running').open('x') as f:f.write('PDAC three-cohort descriptive source analysis')
 try:
  public=run/'public';public.mkdir(exist_ok=False);panel=json.loads((run/'panel.json').read_text());genes=panel['genes'];assert len(genes)==len(set(genes))==173
  outputs=[];cellcounts=[];checks=[]
  for cohort,reader in [('GSE263733',dense263),('GSE278688',sparse278),('GSE242230',sparse242)]:
   print('START',cohort,flush=True);meta,a,lib,present=reader(genes);out,cells,check=summarize(cohort,meta,a,lib,present,genes,run);outputs.append(out);cellcounts.append(cells);checks.append(check)
   out.to_csv(public/(cohort+'_cell_source.tsv'),sep='\t',index=False,na_rep='NA');print('DONE',cohort,json.dumps(check),flush=True)
  allrows=pd.concat(outputs);allrows.to_csv(public/'cell_source_all.tsv',sep='\t',index=False,na_rep='NA');pd.concat(cellcounts).to_csv(public/'celltype_coverage.tsv',sep='\t',index=False)
  consensus=[]
  for gene in genes:
   r={'gene':gene};tops=[]
   for cohort in ['GSE263733','GSE278688','GSE242230']:
    top=allrows[(allrows.gene==gene)&(allrows.cohort==cohort)&(allrows.annotation_level=='broad')&(allrows.source_rank==1)].celltype.tolist();r[cohort]=';'.join(top) or 'NOT_EVALUABLE';tops.append(top[0] if len(top)==1 else None)
   valid=[t for t in tops if t];agreement=max([valid.count(t) for t in set(valid)] or [0]);r['n_evaluable_cohorts']=len(valid);r['max_top_agreement']=agreement;r['interpretation']='THREE_COHORT_DESCRIPTIVE_CONCORDANCE' if agreement==3 else 'TWO_COHORT_DESCRIPTIVE_CONCORDANCE' if agreement==2 else 'COHORT_DEPENDENT_OR_NOT_EVALUABLE';consensus.append(r)
  pd.DataFrame(consensus).to_csv(public/'cross_cohort_source.tsv',sep='\t',index=False)
  SOURCE_FILES.update([run/'panel.json',Path(__file__)])
  pd.DataFrame([{'path':str(p),'sha256':sha(p)} for p in sorted(SOURCE_FILES)]).to_csv(public/'source_manifest.tsv',sep='\t',index=False)
  save(public/'validation.json',{'status':'PASS','cohorts':checks,'patient_rows_exported':False,'library_denominator':'All supplied genes;same cells','raw_count_nonnegative_and_dimensions_checked':True,'minimum_cells_per_unit_type':20,'minimum_units':3,'minimum_units_detected_for_ranking':3,'no_new_p_q':True,'unverified':['Independent re-clustering and malignant CNV labels','Protected patient identity level cross-study de-duplication','Complete biochemistry mapping for all51 metabolites']})
  save(public/'analysis_spec.json',{'code_commit':args.code_commit,'version':'PDAC_sc_source_three_cohorts_v1','panel_genes':genes,'panel_scope':panel,'seed':'NA;deterministic','normalization':'Patient/sample-by-celltype sum counts;log1p(CPM);equal-unit mean;no imputation or batch correction','ranking':'at least20 cells/unit/type;at least3 units and3 detected units/type;descending mean_unit_log1p_cpm','P_q':'Not calculated;descriptive expression sources','versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__}})
  save(run/'DONE.json',{'scope':'Three-cohort173-gene descriptive expression source;51 mapping incomplete'});(run/'.running').unlink()
 except Exception:
  (run/'FAILED.txt').write_text(traceback.format_exc());raise
if __name__=='__main__':main()
