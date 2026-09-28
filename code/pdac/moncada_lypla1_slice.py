"""One-section descriptive spatial display, with honest image-registration status."""
import argparse,json,hashlib,gzip,shutil,csv
from pathlib import Path
import numpy as np,pandas as pd
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize

COLORS={'Cancer':'#ba304d','Duct Epithelium':'#378cb0','Pancreatic':'#d6aa41','Stroma':'#68a871','Unannotated':'#adadad'}
CN={'Cancer':'癌区','Duct Epithelium':'导管上皮区','Pancreatic':'胰腺组织区','Stroma':'间质区','Unannotated':'未注释'}

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()

def setup_axes(ax,title):
 ax.set_xlim(1,33);ax.set_ylim(35,1);ax.set_aspect('equal');ax.set_title(title,fontsize=12);ax.axis('off')

def draw(ax,spots,values,title,vmax,cmap='magma',region=False):
 setup_axes(ax,title)
 patches=PatchCollection([Circle((x,y),.25) for x,y in spots[['x','y']].to_numpy()],edgecolor='none')
 if region:patches.set_facecolor([COLORS[x] for x in values])
 else:patches.set_cmap(cmap);patches.set_norm(Normalize(0,vmax));patches.set_array(np.asarray(values))
 ax.add_collection(patches);return patches

def main():
 p=argparse.ArgumentParser();p.add_argument('--data',required=True);p.add_argument('--run',required=True);p.add_argument('--image',required=True);p.add_argument('--font');args=p.parse_args()
 root=Path(args.run);assert (root/'.running').exists();pub=root/'public';priv=root/'private';pub.mkdir(exist_ok=True);priv.mkdir(exist_ok=True);d=Path(args.data)
 spec=json.loads((root/'moncada_lypla1_slice_spec.json').read_text());shutil.copy(root/'moncada_lypla1_slice_spec.json',pub/'analysis_spec.json')
 x=pd.read_csv(d/'GSM3036911_PDAC-A-ST1-filtered.txt.gz',sep='\t',index_col=0)
 raw=pd.read_csv(d/'GSM3036911.tsv.gz',sep='\t',index_col=0)
 assert x.columns.is_unique and raw.index.is_unique and raw.columns.is_unique
 duplicate_gene_labels=x.index[x.index.duplicated()].tolist()
 assert all(sum(x.index==g)==1 for g in ['LYPLA1','EPCAM','KRT19','COL1A1','PTPRC'])
 assert x.shape==(19738,428) and (x.to_numpy()>=0).all() and np.isfinite(x.to_numpy()).all()
 assert set(x.columns)<=set(raw.index) and 'LYPLA1' in x.index and 'LYPLA1' in raw.columns
 with gzip.open(d/'GSM3036911_PDAC-A-ST1-filtered.txt.gz','rt') as f:
  reader=csv.reader(f,delimiter='\t');header=next(reader);target=[line for line in reader if line[0]=='LYPLA1'];assert len(target)==1
  np.testing.assert_array_equal(np.array(target[0][1:],dtype=float),x.loc['LYPLA1'].to_numpy());assert header[1:]==list(x.columns)
 card=pd.read_csv(d/'card_count_check.tsv',sep='\t',index_col=0)
 np.testing.assert_array_equal(card.loc[x.columns,'LYPLA1'].to_numpy(),raw.loc[x.columns,'LYPLA1'].to_numpy())
 np.testing.assert_array_equal(card.loc[x.columns,'total'].to_numpy(),raw.loc[x.columns].sum(axis=1).to_numpy())
 version_diff=int((x.loc['LYPLA1']!=raw.loc[x.columns,'LYPLA1']).sum())
 spots=pd.DataFrame([list(map(int,k.split('x'))) for k in x.columns],columns=['x','y'],index=x.columns)
 labels=pd.read_csv(d/'layer_manual_PDAC.tsv',sep='\t',index_col=0);assert labels.index.is_unique and len(labels)==426 and set(labels.index)<set(x.columns)
 np.testing.assert_array_equal(spots.loc[labels.index,['x','y']].to_numpy(),labels[['x','y']].to_numpy())
 spots['Region']=labels.Region.reindex(spots.index).fillna('Unannotated');spots['total_counts']=x.sum(axis=0);assert (spots.total_counts>0).all()
 genes=['LYPLA1','EPCAM','KRT19','COL1A1','PTPRC']
 for g in genes:
  assert g in x.index;spots[g+'_count']=x.loc[g];spots[g+'_log1p']=np.log1p(1e4*x.loc[g]/spots.total_counts)
 spots['LYPLA1_fullgene_log1p']=np.log1p(1e4*raw.loc[x.columns,'LYPLA1']/raw.loc[x.columns].sum(axis=1))
 spots.to_csv(priv/'spot_expression_and_regions.tsv.gz',sep='\t',compression='gzip')
 summary=[]
 for region in COLORS:
  z=spots[spots.Region==region]
  summary.append(dict(cancer='PDAC',cohort='GSE111672',sample='GSM3036911',gene='LYPLA1',region=region,spots=len(z),positive_spots=int((z.LYPLA1_count>0).sum()),positive_percent=100*(z.LYPLA1_count>0).mean(),raw_mean=z.LYPLA1_count.mean(),mean_log1p=z.LYPLA1_log1p.mean(),median_log1p=z.LYPLA1_log1p.median(),full_gene_mean_log1p=z.LYPLA1_fullgene_log1p.mean(),p_value='NA',q_value='NA',reason='One patient; descriptive spatial spots, not independent subjects'))
 pd.DataFrame(summary).to_csv(pub/'region_summary.tsv',sep='\t',index=False)
 if args.font:
  font_manager.fontManager.addfont(args.font);plt.rcParams['font.family']=font_manager.FontProperties(fname=args.font).get_name()
 plt.rcParams.update({'axes.unicode_minus':False,'pdf.fonttype':42})
 Image.MAX_IMAGE_PIXELS=None
 # JPEG draft decoding reduces memory; original source remains unchanged.
 with gzip.open(args.image,'rb') as f:
  im=Image.open(f);original_size=im.size;im.draft('RGB',(3500,3500));im.load();im.thumbnail((2500,2500));he=np.asarray(im.convert('RGB'))
 # Consume to EOF to verify gzip CRC, beyond merely decoding a partial JPEG.
 with gzip.open(args.image,'rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''):pass
 fig,axes=plt.subplots(1,3,figsize=(17,7));axes[0].imshow(he);axes[0].axis('off');axes[0].set_title('同一标本的原始 H&E 切片\n完整视野；未与右侧点图作像素配准',fontsize=12)
 draw(axes[1],spots,spots.Region,'已发表区域再注释\nCARD；426 个有标签点',1,region=True)
 a=draw(axes[2],spots,spots.LYPLA1_log1p,'LYPLA1 实测空间分布\n与中间图使用完全相同的点位坐标',float(spots.LYPLA1_log1p.max()));fig.colorbar(a,ax=axes[2],fraction=.033,pad=.02,label='log1p(每万计数)')
 for region,color in COLORS.items():axes[1].scatter([],[],c=color,label=CN[region],s=35)
 axes[1].legend(loc='lower center',bbox_to_anchor=(.5,-.12),ncol=3,frameon=False,fontsize=9)
 fig.suptitle('PDAC-A ST1 / GSE111672：真实切片、区域标签与 LYPLA1',fontsize=18)
 fig.text(.5,.028,'每个点捕获多个细胞；仅一位患者。中、右图可逐点对应；左图目前仅作组织学参照，不能当作已配准的表达叠图。',ha='center',fontsize=11)
 fig.subplots_adjust(left=.01,right=.98,top=.84,bottom=.16,wspace=.13)
 for ext in ['png','pdf']:fig.savefig(pub/('LYPLA1_HE_regions_spatial.'+ext),dpi=190)
 plt.close(fig)
 fig,axes=plt.subplots(2,4,figsize=(17,10))
 draw(axes[0,0],spots,spots.Region,'区域再注释',1,region=True)
 panels=[(spots.LYPLA1_count,'LYPLA1 计数'),(spots.LYPLA1_log1p,'LYPLA1 标准化'),(spots.LYPLA1_fullgene_log1p,'LYPLA1 较早计数版本敏感性')]
 panels.extend((spots[g+'_log1p'],g+' 参考标记') for g in ['EPCAM','KRT19','COL1A1','PTPRC'])
 for ax,(values,title) in zip(axes.ravel()[1:],panels):
  a=draw(ax,spots,values,title,float(max(values.max(),1)));fig.colorbar(a,ax=ax,fraction=.035,pad=.01)
 fig.suptitle('同一组空间坐标：LYPLA1、上皮／间质／免疫参考标记',fontsize=17);fig.text(.5,.02,'各面板标尺分别标明；无插值或平滑。参考标记不等于恶性身份；不把空间点数当作患者数。',ha='center',fontsize=11);fig.subplots_adjust(top=.92,bottom=.065,hspace=.14,wspace=.15)
 for ext in ['png','pdf']:fig.savefig(pub/('LYPLA1_marker_reference_maps.'+ext),dpi=190)
 plt.close(fig)
 sources=[]
 for n,url in [('GSM3036911_PDAC-A-ST1-filtered.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3036nnn/GSM3036911/suppl/GSM3036911_PDAC-A-ST1-filtered.txt.gz'),('GSM3036911.tsv.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3036nnn/GSM3036911/suppl/GSM3036911.tsv.gz'),('Figure4A_layer_annote.RData','https://github.com/YingMa0107/CARD/blob/2d64b91abb5cdd0c7f576b1c5d4727c84e7c93a0/data/Figure4A_layer_annote.RData'),('layer_manual_PDAC.tsv','Derived unchanged from CARD RData via base R write.table')]:sources.append(dict(file=n,source_url=url,bytes=(d/n).stat().st_size,sha256=sha(d/n)))
 sources.append(dict(file=Path(args.image).name,source_url='https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3036nnn/GSM3036911/suppl/GSM3036911_PDAC-A-ST1-HE.jpg.gz',bytes=Path(args.image).stat().st_size,sha256=sha(args.image)))
 pd.DataFrame(sources).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
 val=dict(status='PARTIAL',completed='Actual same-section H&E and measured spatial maps with published regional labels',not_completed='Exact pixel-level H&E-expression registration; multi-patient spatial validation',spots=428,labelled_spots=426,unlabelled_spots=2,patients=1,
  unique_spot_and_plotted_gene_keys=True,unrelated_duplicate_gene_labels_preserved_as_separate_rows=duplicate_gene_labels,label_xy_exact_match=True,independent_filtered_target_read=True,earlier_GEO_counts_and_totals_match_CARD=True,two_GEO_exports_target_discrepant_spots=version_diff,filtered_target_sum=int(x.loc['LYPLA1'].sum()),earlier_target_sum=int(raw.loc[x.columns,'LYPLA1'].sum()),original_he_size=original_size,he_gzip_crc_pass=True,interpolation=False,malignant_single_cell_claim=False)
 (pub/'validation.json').write_text(json.dumps(val,indent=2));(root/'.running').unlink();print(pd.DataFrame(summary).to_string(index=False),flush=True)
if __name__=='__main__':main()
