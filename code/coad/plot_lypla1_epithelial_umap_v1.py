"""Server-only epithelial LYPLA1 UMAP with previously admitted author CNA labels."""
import argparse,json,hashlib,platform
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from run_uckl1_cna_v1 import build_labels
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
PREV=ROOT/'results/collaborative/COAD/B/20260928T031628Z_lypla1_umap_v1'
ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);ap.add_argument('--commit',required=True);a=ap.parse_args()
out=a.out;assert out.parent==PREV.parent and (out/'.running').is_dir();pub=out/'public';pub.mkdir()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
j,paths,coverage,audit=build_labels()
f=pd.read_csv(PREV/'private_lypla1_overlay.tsv',sep='\t',index_col=0,keep_default_na=False)
assert f.index.is_unique and j.cell_id.is_unique and set(j.cell_id)==set(f.index)
j=j.set_index('cell_id');assert (f.case_id==j.loc[f.index,'case_id']).all()
f['group']=j.loc[f.index,'group'];epi=f.loc[f.main_cell_type=='Epithelial'].copy()
assert len(epi)==31247 and np.isfinite(epi[['UMAP1','UMAP2','expression']]).all().all()
groups=['tumor_CNA','tumor_CNN','normal_reference']
expected={'tumor_CNA':4477,'tumor_CNN':7885,'normal_reference':16138}
for g,n in expected.items():assert (epi.group==g).sum()==n
rows=epi.groupby('group').agg(cells=('expression','size'),donors=('case_id','nunique'),detected_cells=('detected','sum'),mean_log1pCP10K=('expression','mean')).reset_index()
rows.to_csv(pub/'coverage.tsv',sep='\t',index=False)
epi.to_csv(out/'private_epithelial_overlay.tsv',sep='\t')
names={'tumor_CNA':'Tumor CNA','tumor_CNN':'Tumor CNN','normal_reference':'Normal reference','excluded_donor_identity_conflict':'Excluded: donor ID conflict','not_selected':'Other / not admitted'}
colors={'tumor_CNA':'#d1495b','tumor_CNN':'#4477aa','normal_reference':'#44aa99','excluded_donor_identity_conflict':'#aaaaaa','not_selected':'#dddddd'}
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
xlim=(epi.UMAP1.min()-.7,epi.UMAP1.max()+.7);ylim=(epi.UMAP2.min()-.7,epi.UMAP2.max()+.7)
vmax=float(epi.loc[epi.group.isin(groups),'expression'].max());norm=Normalize(0,vmax)
def clean(ax,title):ax.set(xlim=xlim,ylim=ylim,xlabel='UMAP 1',ylabel='UMAP 2',xticks=[],yticks=[],title=title);ax.set_aspect('equal')
def bg(ax):ax.scatter(epi.UMAP1,epi.UMAP2,s=1,c='#eeeeee',linewidths=0,rasterized=True)
def expr(ax,z,title):
    bg(ax);zero=z.loc[~z.detected];ax.scatter(zero.UMAP1,zero.UMAP2,s=2,c='#bdbdbd',linewidths=0,rasterized=True)
    p=z.loc[z.detected].sort_values('expression',kind='stable');ax.scatter(p.UMAP1,p.UMAP2,c=p.expression,s=3,cmap='magma',norm=norm,linewidths=0,rasterized=True);clean(ax,title)
def save(fig,name):
    fig.savefig(pub/(name+'.png'),dpi=250,bbox_inches='tight',facecolor='white');p=pub/(name+'.svg');fig.savefig(p,dpi=180,bbox_inches='tight',facecolor='white');plt.close(fig)
    p.write_text('\n'.join(x.rstrip() for x in p.read_text().splitlines())+'\n')
fig,axes=plt.subplots(1,2,figsize=(13,6.5));z=epi.sample(frac=1,random_state=20260928)
axes[0].scatter(z.UMAP1,z.UMAP2,c=z.group.map(colors),s=2,linewidths=0,rasterized=True)
for g in groups+['excluded_donor_identity_conflict','not_selected']:
    n=int((epi.group==g).sum())
    if n:axes[0].scatter([],[],c=colors[g],s=20,label='%s (n=%s)'%(names[g],format(n,',')))
clean(axes[0],'Author CNA groups | epithelial cells');axes[0].legend(loc='lower center',bbox_to_anchor=(.5,-.27),frameon=False,fontsize=8,ncol=2)
expr(axes[1],epi.loc[epi.group.isin(groups)],'LYPLA1 | admitted epithelial cells')
fig.subplots_adjust(left=.05,right=.88,bottom=.25,top=.85,wspace=.17);cax=fig.add_axes([.91,.34,.015,.4]);fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='magma'),cax=cax,label='LYPLA1 log1p(CP10K)')
fig.suptitle('Uhlitz / GSE166555 | epithelial LYPLA1 UMAP',fontsize=16)
fig.text(.05,.025,'Frozen whole-cohort UMAP, cropped to epithelium; no new embedding or batch correction.\nCNN is not proven benign. Excluded / unassigned cells remain pale background in expression views.',fontsize=9);save(fig,'LYPLA1_epithelial_UMAP')
fig,axes=plt.subplots(1,3,figsize=(16,5.7))
for ax,g in zip(axes,groups):expr(ax,epi.loc[epi.group==g],'%s | n=%s'%(names[g],format(expected[g],',')))
fig.subplots_adjust(left=.04,right=.9,bottom=.2,top=.86,wspace=.15);cax=fig.add_axes([.925,.29,.013,.45]);fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='magma'),cax=cax,label='LYPLA1 log1p(CP10K)')
fig.suptitle('LYPLA1 in epithelial cells | same coordinates and color scale',fontsize=16)
fig.text(.04,.035,'Gray: zero counts in the selected group. Pale background: other epithelial cells.\nHigh expression drawn last; pooled visualization only, no new P/q. CNA/CNN are author RNA-inferred copy-number labels.',fontsize=9);save(fig,'LYPLA1_epithelial_UMAP_groups')
sources=[PREV/'private_lypla1_overlay.tsv',paths['metadata'],paths['calls'],Path(__file__),Path(__file__).with_name('run_uckl1_cna_v1.py')]
pd.DataFrame([dict(source_id=p.name,server_path=str(p),sha256=sha(p)) for p in sources]).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
spec=dict(analysis_version='COAD_LYPLA1_epithelial_UMAP_v1',code_lock_commit=a.commit,coordinate_run='20260928T025513Z_uhlitz_umap_v1',expression_run=PREV.name,coordinates_recomputed=False,expression_scale='log1p(CP10K)',vmin=0,vmax=vmax,color_clipping=False,group_rule='reuse run_uckl1_cna_v1.build_labels;exclude whole conflicting donor;normal reference excludes TC labels',excluded_display='annotation explicit;expression pale background only',new_tests=False,software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__))
(pub/'analysis_spec.json').write_text(json.dumps(spec,indent=2)+'\n');(pub/'label_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
v=dict(status='PASS',epithelial_cells=len(epi),admitted_cells=sum(expected.values()),group_counts_match_previous=True,cell_ids_matched=True,donor_ids_matched=True,coordinates_reused=True,new_P_q=False,individual_records_exported=False)
(pub/'validation.json').write_text(json.dumps(v,indent=2)+'\n');(out/'DONE').write_text('DONE\n');(out/'.running').rmdir();print(json.dumps(v),flush=True)
