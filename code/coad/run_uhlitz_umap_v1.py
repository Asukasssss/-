"""Server-only de novo UMAP; retain original author annotations, no new clustering."""
import argparse,gzip,hashlib,io,json,platform,tarfile,time
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
DATA=ROOT/'data/candidates/coad_uhlitz_20260921'
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--commit',required=True);ap.add_argument('--render-only',action='store_true');a=ap.parse_args();out=a.out;pub=out/'public'
assert out.parent==ROOT/'results/collaborative/COAD/B' and len(a.commit)==40
if not a.render_only:assert (out/'.running').is_dir();pub.mkdir()
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def log(x):print(time.strftime('%H:%M:%S'),x,flush=True)
cols=['cell','cell_id','case_id','sample_id','sample_origin','main_cell_type','cell_type_epi_simple','cell_type_imm_simple','cell_type_str_simple','nCount_RNA']
m=pd.read_csv(DATA/'metadata.tsv.gz',sep='\t',usecols=cols,dtype=str,keep_default_na=False)
assert len(m)==68702 and m.cell.is_unique and (m.cell==m.cell_id).all()
m=m.set_index('cell');m['nCount_RNA']=pd.to_numeric(m.nCount_RNA)
assert set(m.sample_origin)=={'Tumor','Normal'}
spec=dict(analysis_version='COAD_Uhlitz_UMAP_v1',code_lock_commit=a.commit,cohort='Uhlitz_GSE166555',coordinate_provenance='New project UMAP; NOT author UMAP coordinates',
    input_scope='All author-retained 68702 cells, 12 donors, tumor and normal; no additional cell filtering',
    annotation='Original author labels; epithelial kept broad in overview; fibroblast subtypes collapsed explicitly; no new clustering or marker inference',
    preprocessing=dict(gene_filter_min_cells=3,normalize_total=10000,log1p=True,HVG_flavor='seurat',HVG_n=2000,scale_zero_center=True,scale_max_value=10),
    embedding=dict(PCA_components=50,PCA_solver='randomized',neighbors=15,metric='euclidean',min_dist=.5,spread=1.,random_state=20260928),
    batch_correction='NONE; donor and tissue views supplied for interpretation',scope_boundary='Visualization only;no DE/P/q;no new malignancy calls;source and per-cell coordinates stay server165')
if not a.render_only:
    hashes={'counts.tar':sha(DATA/'counts.tar'),'metadata.tsv.gz':sha(DATA/'metadata.tsv.gz')}
    assert hashes['counts.tar']=='c5f7c1927be41f27778f5550de46c25b55bdabd60066b87ba4a6bd0d2c6e6b76'
    assert hashes['metadata.tsv.gz']=='4aba042a64a71507db4fad088c8a03510a50a4e9b4202b9bab8d6b07e15303d5'
    blocks=[];cells=[];genes=None;audits=[]
    with tarfile.open(DATA/'counts.tar') as tar:
        entries=[x for x in tar.getmembers() if x.isfile()]
        for i,entry in enumerate(entries):
            stream=io.TextIOWrapper(gzip.GzipFile(fileobj=tar.extractfile(entry)))
            header=stream.readline().rstrip('\r\n').split('\t');assert header[0]=='gene'
            ids=header[1:];assert len(set(ids))==len(ids) and set(ids)<=set(m.index)
            names=[];values=[];indices=[];indptr=[0]
            for line in stream:
                name,text=line.rstrip('\r\n').split('\t',1);v=np.fromstring(text,sep='\t',dtype=np.float64)
                assert len(v)==len(ids) and np.isfinite(v).all() and (v>=0).all() and np.array_equal(v,np.floor(v))
                ix=np.flatnonzero(v);indices.append(ix.astype(np.int32));values.append(v[ix].astype(np.float32));indptr.append(indptr[-1]+len(ix));names.append(name)
            stream.close()
            assert len(names)==len(set(names))
            if genes is None:genes=names
            else:assert genes==names
            block=sp.csr_matrix((np.concatenate(values),np.concatenate(indices),np.asarray(indptr,dtype=np.int64)),shape=(len(names),len(ids))).T.tocsr()
            assert np.array_equal(np.asarray(block.sum(axis=1)).ravel(),m.loc[ids,'nCount_RNA'].to_numpy())
            blocks.append(block);cells.extend(ids);audits.append(dict(block_number=i+1,cells=len(ids),genes=len(names),UMI_matches_metadata=True))
            log('counts block %d/%d admitted'%(i+1,len(entries)))
    assert len(cells)==len(set(cells))==len(m) and set(cells)==set(m.index)
    X=sp.vstack(blocks,format='csr');del blocks
    ad=sc.AnnData(X,obs=m.loc[cells].copy(),var=pd.DataFrame(index=genes));del X
    raw_uckl1=ad[:,['UCKL1']].X.toarray().ravel()
    ad.obs['UCKL1_log1pCP10K']=np.log1p(raw_uckl1/ad.obs.nCount_RNA.to_numpy()*10000)
    sc.settings.n_jobs=4;log('normalize and select HVGs')
    sc.pp.filter_genes(ad,min_cells=3);n_genes_filtered=ad.n_vars
    sc.pp.normalize_total(ad,target_sum=10000);sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad,flavor='seurat',n_top_genes=2000)
    pd.DataFrame(dict(gene=ad.var_names,highly_variable=ad.var.highly_variable.to_numpy())).to_csv(pub/'gene_selection.tsv',sep='\t',index=False)
    ad=ad[:,ad.var.highly_variable].copy();sc.pp.scale(ad,zero_center=True,max_value=10)
    log('PCA');sc.tl.pca(ad,n_comps=50,svd_solver='randomized',random_state=20260928)
    log('neighbors');sc.pp.neighbors(ad,n_neighbors=15,n_pcs=50,metric='euclidean',random_state=20260928)
    log('UMAP');sc.tl.umap(ad,min_dist=.5,spread=1.,random_state=20260928)
    xy=ad.obsm['X_umap'];assert xy.shape==(len(m),2) and np.isfinite(xy).all()
    frame=ad.obs.copy();frame['UMAP1']=xy[:,0];frame['UMAP2']=xy[:,1]
    frame.to_csv(out/'private_cell_embedding.tsv',sep='\t',index=True)
    np.savez_compressed(out/'private_pca.npz',X_pca=ad.obsm['X_pca'],variance_ratio=ad.uns['pca']['variance_ratio'])
    pd.DataFrame(audits).to_csv(pub/'matrix_validation.tsv',sep='\t',index=False)
    source=[]
    for f,url in [('counts.tar','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166555/suppl/GSE166555_RAW.tar'),('metadata.tsv.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166555/suppl/GSE166555_meta_data.tsv.gz')]:source.append(dict(source_id=f,server_path=str(DATA/f),url=url,sha256=hashes[f]))
    pd.DataFrame(source).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
    spec.update(genes_source=len(genes),genes_after_min_cells=n_genes_filtered,actual_HVGs=ad.n_vars,PCA_variance_ratio_sum=float(ad.uns['pca']['variance_ratio'].sum()))
    (pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2)+'\n')
else:
    frame=pd.read_csv(out/'private_cell_embedding.tsv',sep='\t',index_col=0,keep_default_na=False)
    spec=json.loads((pub/'analysis_spec.json').read_text())
frame['display_type']=''
labelmap={'Epithelial':'Epithelial cells','B cells':'B cells','GC B cells':'B cells','Plasma K':'Plasma cells','Plasma L':'Plasma cells','Tconv':'Conventional T cells','CD8+ T':'CD8+ T cells','Treg':'Regulatory T cells','Monocyte':'Monocytes','Mast':'Mast cells','Endothelial':'Endothelial cells','CB FBs':'Fibroblasts','MyoFBs':'Fibroblasts','CAFs':'Fibroblasts','CCL FBs':'Fibroblasts','UC FBs':'Fibroblasts','Pericytes':'Pericytes','Glia':'Glial cells'}
frame['original_label']=np.where(frame.main_cell_type=='Epithelial','Epithelial',np.where(frame.main_cell_type=='Immune',frame.cell_type_imm_simple,frame.cell_type_str_simple))
assert set(frame.original_label)<=set(labelmap)
frame['display_type']=frame.original_label.map(labelmap)
pd.DataFrame([dict(original_label=k,display_label=v,cells=int((frame.original_label==k).sum())) for k,v in labelmap.items()]).to_csv(pub/'annotation_mapping.tsv',sep='\t',index=False)
coverage=frame.groupby(['display_type','sample_origin']).agg(cells=('case_id','size'),donors=('case_id','nunique')).reset_index();coverage.to_csv(pub/'coverage.tsv',sep='\t',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
order=['Epithelial cells','Conventional T cells','CD8+ T cells','Regulatory T cells','B cells','Plasma cells','Monocytes','Mast cells','Fibroblasts','Endothelial cells','Pericytes','Glial cells']
colors=dict(zip(order,['#d1495b','#4477aa','#66aadd','#228833','#aa4499','#cc99bb','#ee7733','#ccbb44','#44aa99','#332288','#999933','#777777']))
rng=np.random.RandomState(20260928);z=frame.iloc[rng.permutation(len(frame))]
def clean(ax):
    ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2');ax.set_xticks([]);ax.set_yticks([]);ax.set_aspect('equal',adjustable='datalim')
def export(fig,name):
    fig.savefig(pub/(name+'.png'),dpi=240,bbox_inches='tight',facecolor='white');fig.savefig(pub/(name+'.svg'),dpi=180,bbox_inches='tight',facecolor='white');plt.close(fig)
fig,ax=plt.subplots(figsize=(12,8))
ax.scatter(z.UMAP1,z.UMAP2,c=z.display_type.map(colors),s=1.3,alpha=.65,linewidths=0,rasterized=True)
for label in order:ax.scatter([],[],s=32,color=colors[label],label='%s  (n=%s)'%(label,format(int((frame.display_type==label).sum()),',')))
clean(ax);ax.legend(loc='center left',bbox_to_anchor=(1.02,.5),frameon=False,fontsize=10,handletextpad=.5)
ax.set_title('Uhlitz / GSE166555 | author cell-type annotations',loc='left',fontsize=16,pad=18)
fig.text(.125,.025,'68,702 cells | 12 donors | tumor + normal tissue\nUMAP recomputed from counts; author labels retained; no batch correction.',fontsize=10,color='#444444')
fig.subplots_adjust(bottom=.13,right=.73);export(fig,'UMAP_celltypes')
fig,axes=plt.subplots(1,2,figsize=(13,6))
for ax,col,palette,title in [(axes[0],'main_cell_type',{'Epithelial':'#d1495b','Immune':'#4477aa','Stromal':'#44aa99'},'Major lineages'),(axes[1],'sample_origin',{'Tumor':'#d1495b','Normal':'#4477aa'},'Tissue of origin')]:
    ax.scatter(z.UMAP1,z.UMAP2,c=z[col].map(palette),s=1,alpha=.6,linewidths=0,rasterized=True)
    for label,color in palette.items():ax.scatter([],[],s=24,c=color,label=label)
    clean(ax);ax.legend(frameon=False,loc='upper right');ax.set_title(title)
fig.suptitle('Uhlitz / GSE166555 | same recomputed UMAP coordinates',fontsize=15);fig.tight_layout();export(fig,'UMAP_lineage_tissue')
# Donor distribution is an interpretation control, not evidence of batch removal.
fig,ax=plt.subplots(figsize=(10,7));donors=sorted(frame.case_id.unique());pal=dict(zip(donors,plt.cm.tab20(np.linspace(0,1,len(donors)))))
ax.scatter(z.UMAP1,z.UMAP2,c=[pal[x] for x in z.case_id],s=1,alpha=.65,linewidths=0,rasterized=True)
for i,label in enumerate(donors):ax.scatter([],[],s=25,color=pal[label],label='Donor %02d'%(i+1))
clean(ax);ax.legend(loc='center left',bbox_to_anchor=(1,.5),frameon=False);ax.set_title('Donor distribution | no batch correction');fig.tight_layout();export(fig,'UMAP_donors')
validation=dict(status='PASS',cells=len(frame),donors=int(frame.case_id.nunique()),tumor_cells=int((frame.sample_origin=='Tumor').sum()),normal_cells=int((frame.sample_origin=='Normal').sum()),display_categories=len(order),coordinate_origin='PROJECT_RECOMPUTED',author_coordinates_available=False,coordinates_finite=bool(np.isfinite(frame[['UMAP1','UMAP2']]).all().all()),cell_ids_unique=bool(frame.index.is_unique),annotation_complete=bool(frame.display_type.ne('').all()),sum_coverage=int(coverage.cells.sum()),no_new_clustering=True,no_new_P_q=True,individual_records_exported=False,software=dict(python=platform.python_version(),scanpy=sc.__version__,numpy=np.__version__,matplotlib=matplotlib.__version__))
assert validation['sum_coverage']==68702 and validation['coordinates_finite'] and validation['annotation_complete']
(pub/'validation.json').write_text(json.dumps(validation,indent=2)+'\n');log(json.dumps(validation))
if not a.render_only:(out/'DONE').write_text('DONE\n');(out/'.running').rmdir()
