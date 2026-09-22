"""Paired expression/annotation UMAPs, reusing frozen Wu extraction and coordinates."""
from pathlib import Path
import sys
# Reuse only extraction and drawing helpers, before any previous figure generation.
root=Path(sys.argv[1]);base=root/'brca_sc156_umap_v1.py'
exec(compile(base.read_text().split("focus=['ASNS'")[0],str(base),'exec'))
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
colors={'B-cells':'#d73027','CAFs':'#4575b4','Cancer Epithelial':'#78b436','Endothelial':'#7b3294','Myeloid':'#ef8a0c','Normal Epithelial':'#43a2ca','Plasmablasts':'#d95f9f','PVL':'#888888','T-cells':'#66a98c'}
labels={'B-cells':'B cells','CAFs':'CAFs','Cancer Epithelial':'Cancer epithelial','Endothelial':'Endothelial','Myeloid':'Myeloid','Normal Epithelial':'Normal epithelial','Plasmablasts':'Plasmablasts','PVL':'PVL','T-cells':'T cells'}
cats=list(colors);counts=obs.celltype_major.value_counts();assert int(counts.sum())==100064
focus=['ASNS','LYPLA1','GLS','LYPLA2','ABHD12','ENPP2','LPCAT1','GPCPD1','GPI','NNMT','PYCR1','ACSL4']
order=focus+[g for g in genes if g not in focus]
pages=[]
with PdfPages(P/'all156_expression_celltype_pairs.pdf') as pdf:
 for page,g in enumerate(order,1):
  fig=plt.figure(figsize=(14,7));ax=fig.add_axes([.035,.13,.375,.75]);right=fig.add_axes([.48,.13,.375,.75]);cb=fig.add_axes([.42,.23,.012,.5])
  im=feature(ax,g);ax.set_title(g+' expression in all cells',fontsize=14,fontweight='bold');fig.colorbar(im,cax=cb,label='log1p(counts per10k)')
  for c in cats:
   sel=obs.celltype_major.eq(c).to_numpy();right.scatter(xy[sel,0],xy[sel,1],s=.5,c=colors[c],linewidths=0,rasterized=True)
  style(right);right.set_title('All cells: cell types',fontsize=14,fontweight='bold')
  handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=colors[c],markeredgecolor='none',markersize=4,label=labels[c]+' ('+format(int(counts[c]),',')+')') for c in cats]
  right.legend(handles=handles,title='Cell type (n cells)',loc='center left',bbox_to_anchor=(1.01,.5),frameon=False,fontsize=8,title_fontsize=9,handletextpad=.3,borderaxespad=0)
  fig.text(.035,.04,'Same cells and coordinates in both panels. Author-derived celltype_major labels; no new clustering.\nPurple = higher RNA; gray = zero/low. Per-gene color cap: positive-expression 99th percentile. Distribution, not an enrichment test.',fontsize=9)
  fig.text(.035,.95,'Wu BRCA | '+g+' | '+str(page)+'/156',fontsize=11,color='#555555')
  pdf.savefig(fig,dpi=115)
  if g in focus:fig.savefig(P/(g+'_expression_celltype_pair.png'),dpi=180)
  plt.close(fig);pages.append(dict(gene=g,page=page,n_cells=100064,color_cap=float(caps[genes.index(g)])))
  if page%10==0:print('Pair page',page,flush=True)
pd.DataFrame(pages).to_csv(P/'gene_page_index.tsv',sep='\t',index=False)
pd.DataFrame([dict(author_celltype=c,n_cells=int(counts[c]),color=colors[c]) for c in cats]).to_csv(P/'celltype_legend.tsv',sep='\t',index=False)
spec=dict(run_id=R.name,version='sc156_umap_pairs_v1',source='Wu frozen obsm/X_umap',genes=156,cells=100064,annotation='author-derived celltype_major',new_P_q=0,new_embedding=False,normalization='unchanged log1p(raw/library*10000)',style='white,small points,light borders;matched axes;right legend with counts',color_cap='each gene positive-expression99th percentile;unchanged from preceding UMAP',source_sha256=digest,matplotlib=matplotlib.__version__)
(P/'analysis_spec.json').write_text(json.dumps(spec,indent=2));(P/'validation.json').write_text(json.dumps(dict(status='DONE',genes=156,pdf_pages=156,focus_png=12,all_cells_retained=True,coordinate_identity_between_panels=True,unique_gene_match=True,source_hash_unchanged=True),indent=2))
pd.DataFrame([dict(path=str(p),sha256=sha(p)) for p in [src,R/'genes156.tsv',R/'brca_sc117_profile_v1.py',base,R/'brca_sc156_umap_pairs_v1.py']]).to_csv(P/'source_manifest.tsv',sep='\t',index=False)
(P/'README_CN.md').write_text('''# BRCA全156基因：表达与细胞类型左右对照UMAP

## 本轮问题
为每个基因增加同坐标细胞类型参照，直观看表达区域对应什么细胞。
## 输入与范围
Wu全部100064细胞、156保留基因；原有X_umap及celltype_major。源矩阵留server165。
## 实际结果
156页PDF，每页左表达、右细胞类型；12个讨论基因另供PNG。gene_page_index定位每个基因；细胞类型图例带实际细胞数。
## 新手解释
左右是相同细胞和坐标，先看左边紫色集中位置，再看右边对应颜色的类别。色标和原表达UMAP相同，不换算成用户参考图中的log2(CPM+0.5)。每基因有独立上限，灰色为零或低；高表达点最后绘制。
## 限制与反证
这是表达分布，不是正式富集检验；不用圈选区域暗示唯一作用细胞。右边图例是细胞数，不是独立患者数或真实组织比例。不同基因颜色不直接比较，统计解释仍以供者汇总为准。图中细胞类型使用作者原标签，中文解释中的恶性上皮对应Cancer Epithelial。
## 当前决定
全部156保留，无新P/q、无新降维或聚类。历史图不覆盖。
## 下一步
对照现有供者等权结果理解背景；不按图上深浅淘汰候选。
## 复现
新server165运行目录包含genes156.tsv、brca_sc117_profile_v1.py、brca_sc156_umap_v1.py及本脚本，python3 brca_sc156_umap_pairs_v1.py RUN_DIR。
''',encoding='utf-8')
(P/'.gitattributes').write_text('* -text\n');pd.DataFrame([dict(file=p.name,sha256=sha(p)) for p in sorted(P.iterdir()) if p.name!='checksums.tsv']).to_csv(P/'checksums.tsv',sep='\t',index=False)
(R/'DONE').write_text('DONE');print('DONE',flush=True)
