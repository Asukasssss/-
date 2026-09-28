"""Server-only LYPLA1 expression overlay on frozen Uhlitz UMAP coordinates."""
import argparse,gzip,hashlib,io,json,platform,tarfile
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
DATA=ROOT/'data/candidates/coad_uhlitz_20260921'
OLD=ROOT/'results/collaborative/COAD/B/20260928T025513Z_uhlitz_umap_v1'
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--commit',required=True);ap.add_argument('--render-only',action='store_true');a=ap.parse_args()
out=a.out;assert out.parent==OLD.parent and (a.render_only or (out/'.running').is_dir())
pub=out/'public'
if not a.render_only:pub.mkdir()
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
if not a.render_only:
    manifest=[]
    for p in [DATA/'counts.tar',DATA/'metadata.tsv.gz',OLD/'private_cell_embedding.tsv']:
        manifest.append(dict(source_id=p.name,server_path=str(p),sha256=sha(p)))
    assert manifest[0]['sha256']=='c5f7c1927be41f27778f5550de46c25b55bdabd60066b87ba4a6bd0d2c6e6b76'
    assert manifest[1]['sha256']=='4aba042a64a71507db4fad088c8a03510a50a4e9b4202b9bab8d6b07e15303d5'
    frame=pd.read_csv(OLD/'private_cell_embedding.tsv',sep='\t',index_col=0,keep_default_na=False)
    assert len(frame)==68702 and frame.index.is_unique and (frame.nCount_RNA>0).all()
    values={};blocks=0
    with tarfile.open(DATA/'counts.tar') as tar:
        for entry in tar.getmembers():
            if not entry.isfile():continue
            with io.TextIOWrapper(gzip.GzipFile(fileobj=tar.extractfile(entry))) as stream:
                header=stream.readline().rstrip('\r\n').split('\t');assert header[0]=='gene';ids=header[1:];found=0
                for line in stream:
                    name,text=line.rstrip('\r\n').split('\t',1)
                    if name!='LYPLA1':continue
                    v=np.fromstring(text,sep='\t',dtype=float)
                    assert len(v)==len(ids) and np.isfinite(v).all() and (v>=0).all() and np.array_equal(v,np.floor(v))
                    for cell,x in zip(ids,v):
                        assert cell not in values
                        values[cell]=x
                    found+=1
                assert found==1
            blocks+=1;print('gene extracted block',blocks,flush=True)
    assert blocks==25 and set(values)==set(frame.index)
    frame['LYPLA1_counts']=[values[x] for x in frame.index]
    frame['expression']=np.log1p(frame.LYPLA1_counts.to_numpy()/frame.nCount_RNA.to_numpy()*10000)
    frame['detected']=frame.LYPLA1_counts>0
    frame.to_csv(out/'private_lypla1_overlay.tsv',sep='\t',index=True)
else:
    frame=pd.read_csv(out/'private_lypla1_overlay.tsv',sep='\t',index_col=0,keep_default_na=False)
    manifest=pd.read_csv(pub/'source_manifest.tsv',sep='\t').to_dict('records')
frame['original_label']=np.where(frame.main_cell_type=='Epithelial','Epithelial',np.where(frame.main_cell_type=='Immune',frame.cell_type_imm_simple,frame.cell_type_str_simple))
mapping=pd.read_csv(OLD/'public/annotation_mapping.tsv',sep='\t').set_index('original_label').display_label
frame['display_type']=frame.original_label.map(mapping);assert frame.display_type.notna().all()
summary=frame.groupby(['sample_origin','display_type']).agg(cells=('expression','size'),donors=('case_id','nunique'),detected_cells=('detected','sum'),mean_log1pCP10K=('expression','mean'),median_log1pCP10K=('expression','median')).reset_index()
summary['detected_fraction']=summary.detected_cells/summary.cells
summary.to_csv(pub/'expression_coverage.tsv',sep='\t',index=False)
pd.DataFrame(manifest).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
vmax=float(frame.expression.max());assert vmax>0
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
xlim=(frame.UMAP1.min()-1,frame.UMAP1.max()+1);ylim=(frame.UMAP2.min()-1,frame.UMAP2.max()+1)
def panel(ax,f,title):
    ax.scatter(frame.UMAP1,frame.UMAP2,s=.7,color='#eeeeee',linewidths=0,rasterized=True)
    zero=f.loc[~f.detected];ax.scatter(zero.UMAP1,zero.UMAP2,s=1,color='#c9c9c9',linewidths=0,rasterized=True)
    z=f.loc[f.detected].sort_values('expression',kind='stable')
    ax.scatter(z.UMAP1,z.UMAP2,c=z.expression,s=2,cmap='magma',norm=Normalize(0,vmax),linewidths=0,rasterized=True)
    ax.set(xlim=xlim,ylim=ylim,xlabel='UMAP 1',ylabel='UMAP 2',xticks=[],yticks=[],title=title);ax.set_aspect('equal')
def save(fig,name):
    fig.savefig(pub/(name+'.png'),dpi=250,bbox_inches='tight',facecolor='white')
    p=pub/(name+'.svg');fig.savefig(p,dpi=180,bbox_inches='tight',facecolor='white');plt.close(fig)
    p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n')
fig,ax=plt.subplots(figsize=(9,8));panel(ax,frame,'LYPLA1 | Uhlitz / GSE166555 | all cells')
fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,vmax),cmap='magma'),ax=ax,shrink=.65,label='LYPLA1 log1p(counts / total UMI x 10,000)')
fig.text(.12,.025,'68,702 cells | 12 donors | tumor + normal\nGray: zero detected counts; high expression drawn last.\nSame project UMAP as annotation figure; no batch correction.',fontsize=9)
fig.subplots_adjust(bottom=.16);save(fig,'LYPLA1_UMAP')
fig,axes=plt.subplots(1,2,figsize=(14,6.5))
for ax,tissue in zip(axes,['Tumor','Normal']):
    f=frame.loc[frame.sample_origin==tissue];panel(ax,f,'%s tissue | %s cells'%(tissue,format(len(f),',')))
fig.subplots_adjust(bottom=.13,top=.9,right=.86,wspace=.2)
cax=fig.add_axes([.89,.25,.018,.5])
fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,vmax),cmap='magma'),cax=cax,label='LYPLA1 log1p(CP10K)')
fig.suptitle('LYPLA1 | same coordinates and color scale',fontsize=15)
fig.text(.12,.025,'Gray: zero counts; light-gray background: other tissue. Tumor origin does not imply malignant identity.',fontsize=9)
save(fig,'LYPLA1_UMAP_tissue')
spec=dict(analysis_version='COAD_LYPLA1_UMAP_v1',code_lock_commit=a.commit,gene='LYPLA1',cohort='Uhlitz_GSE166555',coordinate_run=OLD.name,coordinate_recomputed=False,normalization='log1p(raw gene UMI / author nCount_RNA * 10000)',color_scale=[0,vmax],color_clipping=False,zero_color='gray',draw_order='nonzero ascending, high last',batch_correction='NONE',scope='68702 cells from 12 donors;all author retained cells;not CNA-only subset',statistical_tests='NONE',aggregate_unit='pooled cells, descriptive only;not donor-equally-weighted',software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__))
(pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2)+'\n')
validation=dict(status='PASS',cells=len(frame),donors=int(frame.case_id.nunique()),detected_cells=int(frame.detected.sum()),zero_cells=int((~frame.detected).sum()),all_cells_matched=True,source_hashes_match=True,same_embedding=True,new_P_q=False,individual_records_exported=False)
assert summary.cells.sum()==len(frame) and np.isfinite(frame.expression).all()
(pub/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
if not a.render_only:
    (out/'DONE').write_text('DONE\n');(out/'.running').rmdir()
print(json.dumps(validation),flush=True)
