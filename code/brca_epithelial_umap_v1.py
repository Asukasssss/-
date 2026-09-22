"""New exploratory epithelial embedding, frozen author labels, all156 triptychs."""
from pathlib import Path
import sys
root=Path(sys.argv[1]);base=root/'brca_sc156_umap_v1.py'
exec(compile(base.read_text().split("focus=['ASNS'")[0],str(base),'exec'))
from sklearn.decomposition import PCA
from umap import UMAP
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
import sklearn,umap
keep=obs.celltype_major.isin(['Cancer Epithelial','Normal Epithelial']).to_numpy();assert keep.sum()==28844
epobs=obs.loc[keep].copy();epobs['group']=np.where(epobs.celltype_major.eq('Normal Epithelial'),'Normal epithelial',epobs.celltype_minor)
cats=['Cancer Basal SC','Cancer Cycling','Cancer Her2 SC','Cancer LumA SC','Cancer LumB SC','Normal epithelial']
colors=['#d73027','#4575b4','#78b436','#7b3294','#ef8a0c','#bfc5cc'];counts=epobs.group.value_counts()
assert [int(counts[c]) for c in cats]==[4312,5359,3708,7742,3368,4355]
blocks=[]
with h5py.File(src,'r') as h:
 raw=h['raw/X'];ptr=raw['indptr'][:]
 for start in range(0,len(obs),2000):
  stop=min(start+2000,len(obs));mask=keep[start:stop]
  if not mask.any():continue
  lo,hi=ptr[start],ptr[stop];a=sparse.csr_matrix((raw['data'][lo:hi],raw['indices'][lo:hi],ptr[start:stop+1]-lo),shape=(stop-start,len(var)))[mask].astype(np.float32)
  a=a.multiply(10000/np.asarray(a.sum(axis=1))).tocsr() if False else a.multiply((10000/np.asarray(a.sum(axis=1)).ravel())[:,None]).tocsr()
  a.data=np.log1p(a.data);blocks.append(a)
X=sparse.vstack(blocks,format='csr');del blocks
mu=np.asarray(X.mean(axis=0)).ravel();variance=np.asarray(X.multiply(X).mean(axis=0)).ravel()-mu**2
eligible=np.array([not (str(n).startswith(('MT-','RPL','RPS'))) for n in names]);eligible&=np.asarray((X>0).sum(axis=0)).ravel()>=20
pool=np.flatnonzero(eligible);selected=pool[np.argsort(variance[pool],kind='stable')[-2000:]]
pd.DataFrame(dict(feature_id=var.index[selected],symbol=names[selected],log_expression_variance=variance[selected])).to_csv(P/'embedding_genes.tsv',sep='\t',index=False)
pc=PCA(n_components=30,svd_solver='randomized',random_state=20260922);scores=pc.fit_transform(X[:,selected].toarray());del X
print('PCA done',flush=True)
model=UMAP(n_neighbors=30,min_dist=.3,n_components=2,metric='euclidean',random_state=20260922,n_jobs=1)
xy=model.fit_transform(scores);assert xy.shape==(28844,2) and np.isfinite(xy).all()
private=R/'private';private.mkdir(exist_ok=True);np.save(private/'epithelial_umap.npy',xy)
epobs[['group','donor_id','subtype']].to_csv(private/'epithelial_cell_metadata.tsv',sep='\t')
expr=expr[keep];obs=epobs;caps=np.array([np.quantile(expr[:,i][expr[:,i]>0],.99) if (expr[:,i]>0).any() else 1 for i in range(156)])
limits=[xy[:,0].min()-.4,xy[:,0].max()+.4,xy[:,1].min()-.4,xy[:,1].max()+.4]
focus=['ASNS','LYPLA1','GLS','LYPLA2','ABHD12','ENPP2','LPCAT1','GPCPD1','GPI','NNMT','PYCR1','ACSL4'];order=focus+[g for g in genes if g not in focus]
pages=[]
with PdfPages(P/'all156_epithelial_triptychs.pdf') as pdf:
 for page,g in enumerate(order,1):
  fig=plt.figure(figsize=(16,6));axes=[fig.add_axes(v) for v in [[.025,.15,.25,.70],[.345,.15,.25,.70],[.615,.15,.25,.70]]];cb=fig.add_axes([.282,.25,.01,.5])
  im=feature(axes[0],g);axes[0].set_title(g+' | epithelial expression',fontsize=12);fig.colorbar(im,cax=cb,label='log1p(counts per10k)')
  ax=axes[1];ax.scatter(xy[:,0],xy[:,1],s=.5,c='#e4e4e4',linewidths=0,rasterized=True);sel=obs.group.eq('Cancer Basal SC').to_numpy();ax.scatter(xy[sel,0],xy[sel,1],s=.6,c=colors[0],linewidths=0,rasterized=True);style(ax);ax.set_title('Cancer Basal SC (n=4,312)',fontsize=12)
  ax=axes[2]
  for c,color in zip(cats,colors):
   sel=obs.group.eq(c).to_numpy();ax.scatter(xy[sel,0],xy[sel,1],s=.6,c=color,linewidths=0,rasterized=True)
  style(ax);ax.set_title('Epithelial cells',fontsize=12)
  handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=color,markeredgecolor='none',markersize=4,label=c+' ('+format(int(counts[c]),',')+')') for c,color in zip(cats,colors)]
  ax.legend(handles=handles,title='Author-derived annotation',loc='center left',bbox_to_anchor=(1.01,.5),frameon=False,fontsize=7,title_fontsize=8,handletextpad=.2,borderaxespad=0)
  fig.text(.025,.94,'Wu BRCA | epithelial-only UMAP | '+g+' | '+str(page)+'/156',fontsize=13,fontweight='bold')
  fig.text(.025,.035,'NEW exploratory embedding of all 28,844 epithelial cells; author labels retained, no new clustering.\nSame coordinates in three panels. Color cap per gene: positive-expression 99th percentile. Basal cell-state label is not patient TNBC status.',fontsize=9)
  pdf.savefig(fig,dpi=110)
  if g in focus:fig.savefig(P/(g+'_epithelial_triptych.png'),dpi=180)
  plt.close(fig);pages.append(dict(gene=g,page=page,color_cap=float(caps[genes.index(g)])))
  if page%20==0:print('Page',page,flush=True)
pd.DataFrame(pages).to_csv(P/'gene_page_index.tsv',sep='\t',index=False);pd.DataFrame([dict(celltype=c,n_cells=int(counts[c])) for c in cats]).to_csv(P/'epithelial_group_counts.tsv',sep='\t',index=False)
spec=dict(run_id=R.name,analysis_version='epithelial_umap_v1',source_sha256=digest,scope='all28844 epithelial cells;all156 genes',annotation='frozen celltype_minor;all normal epithelial subtypes pooled',embedding='NEW, not author coordinates',normalization='log1p(raw/all-raw-gene library*10000)',feature_selection='top2000 log-expression variance;detected>=20cells;exclude MT-/RPL/RPS prefix;not dispersion-corrected HVG',PCA='centered unscaled,30 randomized components',pca_variance_explained=float(pc.explained_variance_ratio_.sum()),UMAP=dict(n_neighbors=30,min_dist=.3,metric='euclidean',random_state=20260922,n_jobs=1),new_clustering=False,batch_correction=False,new_P_q=0,software=dict(sklearn=sklearn.__version__,umap=umap.__version__,matplotlib=matplotlib.__version__))
(P/'analysis_spec.json').write_text(json.dumps(spec,indent=2));(P/'validation.json').write_text(json.dumps(dict(status='DONE',all_epithelial_cells_retained=True,cells=28844,genes=156,pages=156,counts_match_reference=True,finite_embedding=True,source_hash_unchanged=digest==sha(src)),indent=2))
pd.DataFrame([dict(path=str(p),sha256=sha(p)) for p in [src,base,R/'brca_sc117_profile_v1.py',R/'genes156.tsv',R/'brca_epithelial_umap_v1.py']]).to_csv(P/'source_manifest.tsv',sep='\t',index=False)
(P/'README_CN.md').write_text('''# 上皮细胞专用UMAP三联图
## 本轮问题
参考用户图，左侧基因表达，中间Cancer Basal SC高亮，右侧全部上皮亚群注释。
## 输入与范围
Wu癌上皮24489+正常上皮4355=28844细胞，156个保留基因。原作者minor标签保留；正常三个亚群合并为Normal epithelial。
## 实际结果
新计算上皮专用UMAP，156页三联图及12个重点PNG。亚群计数与用户示例一致，但坐标不是示例坐标的复刻。
## 新手解释
紫色深表示RNA较高；灰色为零/低；中间红色只是标出基底样癌上皮细胞，不表示该基因只在此表达。Cancer Basal SC、LumA等是细胞注释，不等于患者的临床TNBC、ER+分型。
## 限制与反证
新UMAP仅作探索可视化；使用前2000个log表达方差较大的基因，未做批次校正，形状可能受供者和技术差异影响。未根据岛屿重新聚类或宣称新亚群。图上距离、岛屿或颜色不能证明功能、富集或分型差异。色标仍用log1p(counts per10k)，不是示例的log2(CPM+0.5)；每基因阳性值99分位截断，细胞全部保留。
## 当前决定
作为上皮表达补充图，不替代全细胞图或供者等权分析，不能因此忽略免疫/基质表达。
## 下一步
用于读图，不自动增加机制或亚型统计。
## 复现
server165新目录含基因清单和3个脚本，运行python3 brca_epithelial_umap_v1.py RUN_DIR。随机种子20260922，详细参数与软件版本见analysis_spec；逐细胞坐标和标签仅存服务器private。
''',encoding='utf-8');(P/'.gitattributes').write_text('* -text\n')
pd.DataFrame([dict(file=p.name,sha256=sha(p)) for p in sorted(P.iterdir()) if p.name!='checksums.tsv']).to_csv(P/'checksums.tsv',sep='\t',index=False);(R/'DONE').write_text('DONE');print('DONE',flush=True)
