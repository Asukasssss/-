"""Prepare per-donor inferCNV inputs on server165; all cell data remain private."""
import argparse,gzip,re,json,hashlib,platform
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse,io
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--source-run',type=Path,required=True);ap.add_argument('--code-commit',required=True);a=ap.parse_args()
out=a.out;pub=out/'public/06_EXTERNAL';private=out/'private';src=a.source_run/'source/PRAD24_cellxgene.h5ad'
gtf=Path('/public3/xuzx/Cancer/lung_cancer_pack/raw/bulk_transcriptomics/xena_gdc_tcga/gencode.v36.annotation.gtf.gz')
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def js(p,d):p.write_text(json.dumps(d,indent=2,ensure_ascii=False),encoding='utf-8')
obj=ad.read_h5ad(src,backed='r');o=obj.obs.copy();assert o.index.is_unique
assert o.groupby('sample',observed=True).donor_id.nunique().max()==1
rows=[];header=[]
with gzip.open(gtf,'rt') as f:
 for line in f:
  if line.startswith('#'):header.append(line.strip());continue
  s=line.rstrip().split('\t')
  if s[2]!='gene' or s[0] not in ['chr'+str(i) for i in range(1,23)]:continue
  attrs=dict(re.findall(r'(\w+) "([^"]*)"',s[8]));gid=attrs['gene_id'].split('.')[0]
  rows.append(dict(gene=gid,chr=s[0],start=int(s[3]),end=int(s[4]),symbol=attrs['gene_name']))
assert 'GRCh38' in header[0]
g=pd.DataFrame(rows);dup=g.gene.duplicated(keep=False);ambiguous=g.loc[dup,'gene'].unique().tolist();g=g[~dup].set_index('gene')
ids=pd.Index([str(x).split('.')[0] for x in obj.raw.var.index]);assert ids.is_unique
keep=ids.isin(g.index);ix=np.flatnonzero(keep);genes=ids[keep];g=g.loc[genes];g.index.name='gene';g=g.reset_index();g['chr_num']=g.chr.str[3:].astype(int)
order=g.sort_values(['chr_num','start','end','gene']);ix=ix[order.index.to_numpy()];order=order.reset_index(drop=True)
order[['gene','chr','start','end']].to_csv(out/'source/gene_order.tsv',sep='\t',header=False,index=False)
save(order,pub/'gene_coordinate_dictionary.tsv')
seed=20260925;rng=np.random.default_rng(seed)
plans=[];metadata=[];sels={};excluded=[]
donors=sorted(o.donor_id.astype(str).unique())
for k,donor in enumerate(donors,1):
 alias=f'D{k:02d}';d=o[o.donor_id.astype(str).eq(donor)]
 my=d[d.celltype_major_v2.astype(str).eq('Myeloid')]
 t=d[d.celltype_major_v2.astype(str).eq('T-cells')&~d.celltype_minor_v2.astype(str).isin(['NK_cells','NK_like_T-cells','T-cells_IFN'])]
 if len(my)<10 or len(t)<100:
  excluded.append(dict(donor_alias=alias,reason='fewer_than10_myeloid_or100_Tcells',n_myeloid=len(my),n_T=len(t)));continue
 # Disjoint reference pools; stratification balances tissue contribution when possible.
 poolA=[];poolB=[]
 for tissue,z in t.groupby('type',observed=True):
  names=rng.permutation(z.index.to_numpy());half=len(names)//2
  poolA.extend(names[:min(100,half)]);poolB.extend(names[half:half+min(100,len(names)-half)])
 if min(len(poolA),len(poolB))<50:
  excluded.append(dict(donor_alias=alias,reason='fewer_than50_reference_cells_per_split',n_myeloid=len(my),n_T=len(t)));continue
 b=d[d.celltype_major_v2.astype(str).eq('B-cells')];bs=rng.choice(b.index,min(50,len(b)),replace=False).tolist()
 ep=d[d.celltype_major_v2.astype(str).eq('Epithelial')&d.malignant_anno_merged.astype(str).eq('malignant')&d['type'].astype(str).eq('cancer')]
 eps=rng.choice(ep.index,min(100,len(ep)),replace=False).tolist()
 selected=list(my.index)+poolA+poolB+bs+eps;assert len(set(selected))==len(selected)
 sels[alias]=selected
 for name in selected:
  z=o.loc[name]
  role='Myeloid' if name in my.index else 'T_A' if name in poolA else 'T_B' if name in poolB else 'B_control' if name in bs else 'Epithelial_control'
  group=role if role!='Myeloid' else 'MY_'+str(z['type'])+'_'+str(z.celltype_minor_v2)
  metadata.append(dict(cell=name,donor=donor,donor_alias=alias,role=role,group=group,tissue=str(z['type']),subtype=str(z.celltype_minor_v2)))
 plans.append(dict(donor_alias=alias,n_myeloid=len(my),n_reference_A=len(poolA),n_reference_B=len(poolB),n_B_control=len(bs),n_epithelial_control=len(eps),n_cells=len(selected)))
meta=pd.DataFrame(metadata);save(meta,private/'selected_cell_metadata.tsv');save(pd.DataFrame(plans),private/'donor_plan.tsv');save(pd.DataFrame(excluded),private/'excluded_donors.tsv')
all_selected=set(meta.cell);mask=o.index.isin(all_selected);positions=np.flatnonzero(mask);chunks=[]
for start in range(0,len(o),2000):
 end=min(start+2000,len(o));sel=mask[start:end]
 if not sel.any():continue
 x=obj.raw.X[start:end].tocsr()[sel][:,ix]
 assert np.isfinite(x.data).all() and (x.data>=0).all() and np.equal(x.data,np.round(x.data)).all()
 chunks.append(x)
mat=sparse.vstack(chunks,format='csr');selected_names=o.index[positions];assert len(selected_names)==mat.shape[0]
for alias,names in sels.items():
 sub=private/alias;sub.mkdir(exist_ok=False)
 loc=selected_names.get_indexer(names);assert (loc>=0).all();counts=mat[loc].T.tocsc()
 io.mmwrite(str(sub/'counts.mtx'),counts)
 (sub/'genes.txt').write_text('\n'.join(order.gene)+'\n');(sub/'cells.txt').write_text('\n'.join(names)+'\n')
 z=meta.set_index('cell').loc[names]
 for mode,reference in [('primary','T_A'),('reference_split','T_B')]:
  groups=z.group.copy();groups[z.role.eq(reference)]='REF_T';groups[z.role.eq('T_B' if reference=='T_A' else 'T_A')]='T_holdout'
  pd.DataFrame({'cell':names,'group':groups.to_numpy()}).to_csv(sub/f'annotations_{mode}.tsv',sep='\t',header=False,index=False)
print('INPUTS',len(plans),'donors',len(meta),'cells',len(order),'genes',flush=True)
spec=dict(version='prad_myeloid_infercnv_v1',parent_commit=a.code_commit,input_dataset_version='68b23fda-7191-46a5-8870-819feca3e66e',author_annotations='major_v2 Myeloid plus minor_v2;unchanged',method='infercnv official R package;independent donor runs;two disjoint T-reference splits',target='all myeloid cells in eligible donors, tumor and adjacent',reference='same donor T cells excluding author NK/NK-like/IFN classes;two disjoint pools stratified by tissue;50minimum per pool',controls='other reference split held out;B cells up to50 and author malignant tumor epithelium up to100;controls not used as anchors',minimum_myeloid_per_donor=10,reference_cap_per_tissue_per_pool=100,seed=seed,coordinate_source='GENCODEv36 GRCh38/Ensembl102;exact unique Ensembl stable ID mapping;autosomes only',coordinate_build_scope='reference coordinates chosen explicitly for gene ordering, not new verification of original alignment build',gtf_header=header,raw_gene_count=len(ids),mapped_autosomal_genes=len(order),excluded_ambiguous_gene_ids=ambiguous,eligible_donors=len(plans),excluded_donors=len(excluded),eligible_myeloid_cells=int(sum(p['n_myeloid'] for p in plans)),all_myeloid_cells=int(o.celltype_major_v2.astype(str).eq('Myeloid').sum()),parameters=dict(cutoff=.1,min_cells_per_gene=3,window_length=101,cluster_by_groups=True,denoise=True,HMM=False,num_threads=2),scope='RNA-inferred broad copy-number-like signal, not DNA CNV or malignancy classification;HMM disabled, no discrete CNV calls;disjoint reference sensitivity not independent cohort validation',software_preparation=dict(python=platform.python_version(),anndata=ad.__version__,numpy=np.__version__,pandas=pd.__version__))
js(pub/'analysis_spec.json',spec)
save(pd.DataFrame([dict(source_id=f.name,path_or_url=str(f),sha256=sha(f)) for f in [src,gtf,Path(__file__)]]),pub/'source_manifest.tsv')
js(pub/'input_validation.json',dict(status='DONE',cells_unique=True,genes_unique=True,raw_integer=True,reference_splits_disjoint=True,no_myeloid_reference=True,no_new_annotation=True))
obj.file.close()
