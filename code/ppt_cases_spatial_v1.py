"""Server-only gene extraction on existing spatial sections; no new tests."""
from pathlib import Path
import sys,json,io,zipfile,hashlib
import numpy as np,pandas as pd,h5py
from scipy.sparse import csc_matrix
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from PIL import Image
R=Path(sys.argv[1]);P=R/'public';P.mkdir(exist_ok=True)
B=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716');C=B/'results/collaborative';AS=C/'BRCA/A/20260929T071500Z_four_sc_v3'
for n in ['msyh.ttc','msyhbd.ttc']:font_manager.fontManager.addfont(str(AS/n))
plt.rcParams.update({'font.family':'Microsoft YaHei','pdf.fonttype':42,'axes.unicode_minus':False})
G='#00553B';INK='#183D34';M='#6A7973';LINE='#DDE7E0';CAP=3.
COL={'Tumor':'#C94B54','Cancer':'#C94B54','Benign':'#3276A5','Immune':'#D1AE35','DCIS':'#836BB1','Unannotated':'#BFC6C9','Mixed':'#E4E8E8','Stroma':'#BFC6C9','PIN':'#A98C24','Transition_State':'#258F88','Vessel':'#258F88','Necrosis':'#777777','Fat':'#C4AE94','Nerve':'#749487','Non-neoplastic epithelium':'#3276A5','Tumor-stroma mixed':'#EDAC68','Other annotated tissue':'#A481B1'}
LABEL={'Tumor':'癌区','Cancer':'癌区','Benign':'良性腺体','Immune':'免疫区','DCIS':'原位癌','Unannotated':'未标注区','Mixed':'混合区','Stroma':'间质','PIN':'PIN','Transition_State':'过渡区','Vessel':'血管','Necrosis':'坏死','Fat':'脂肪区','Nerve':'神经区','Non-neoplastic epithelium':'非肿瘤上皮','Tumor-stroma mixed':'癌 / 间质混合','Other annotated tissue':'其他组织'}
sources=[];checks=[];data=[];summ=[]
def source(p):sources.append(dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest() if p.stat().st_size<100000000 else 'large source; inherited manifest'))
def dec(a):return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in a])
def gene10x(path,gene):
 with h5py.File(path) as h:
  m=h['matrix'];x=csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=m['shape'][:]);names=dec(m['features/name'][:]);bc=dec(m['barcodes'][:]);ft=dec(m['features/feature_type'][:]);ix=np.flatnonzero((names==gene)&(ft=='Gene Expression'));assert len(ix)==1
  total=np.asarray(x[ft=='Gene Expression'].sum(0)).ravel();cnt=x[ix[0]].toarray().ravel();vals=np.log1p(np.divide(cnt*10000.,total,out=np.zeros(len(total)),where=total>0))
 return pd.Series(vals,index=bc),pd.Series(total,index=bc)
for s in ['V19T26-012_B1','V19T26-032_B1']:
 r=C/'BRCA/A/20260928T060000Z_cta_lypla1_v1';p=r/(s+'_spots_SERVER_ONLY.tsv.gz');source(p);d=pd.read_csv(p,sep='\t');d=d[d.qc.astype(str).str.lower().isin(['true','1'])]
 folder=B/'data/candidates/BRCA_CTA2025/spaceranger_output'/s/'outs';hp=folder/'filtered_feature_bc_matrix.h5';source(hp);v,total=gene10x(hp,'ASNS');assert set(d.barcode)<=set(v.index);assert np.allclose(total.loc[d.barcode],d.total)
 ip=folder/'spatial/tissue_lowres_image.png';source(ip)
 data.append(dict(cancer='BRCA',gene='ASNS',sample=s,patient=d.patient.iloc[0],img=np.asarray(Image.open(ip)),xy=d[['x','y']].to_numpy(),expr=v.loc[d.barcode].to_numpy(),region=d.region.to_numpy()))
r=C/'PRAD/B/20260928T101700Z_lypla1_erickson_v1';p=r/'private/spot_values.tsv.gz';source(p);d=pd.read_csv(p,sep='\t')
for s in ['H2_5','H2_1']:
 x=d[d.section==s];folder=r/'source/Patient1'/s;hp=folder/'filtered_feature_bc_matrix.h5';source(hp);v,total=gene10x(hp,'SLC6A6');assert np.allclose(total.loc[x.barcode],x.total_umi)
 sf=json.loads((folder/'scalefactors_json.json').read_text());ip=folder/(s+'_tissue_hires_image.png');source(ip)
 data.append(dict(cancer='PRAD',gene='SLC6A6',sample=s,patient='Patient 1',img=np.asarray(Image.open(ip)),xy=x[['pixel_col','pixel_row']].to_numpy()*sf['tissue_hires_scalef'],expr=v.loc[x.barcode].to_numpy(),region=x.region.to_numpy()))
# First two source-design colon/cecum cases, selected without target expression.
for alias,name in [('S1_Cec_Rep1','SN123_A551763_Rep1'),('S2_Col_R_Rep1','SN123_A595688_Rep1')]:
 p=B/'data/candidates/coad_lypla1_pathology_20260928'/(name+'.selected.zip');source(p)
 dp=C/'COAD/B/20260928T094500Z_lypla1_pathology_v1/private'/(alias+'_spots.tsv');source(dp);d=pd.read_csv(dp,sep='\t');d=d[d.admitted.astype(str).str.lower().isin(['true','1'])]
 with zipfile.ZipFile(p) as z:
  def member(suffix):return next(n for n in z.namelist() if n.endswith(suffix) and '__MACOSX' not in n)
  v,total=gene10x(io.BytesIO(z.read(member('filtered_feature_bc_matrix.h5'))),'UCKL1');img=np.asarray(Image.open(io.BytesIO(z.read(member('tissue_lowres_image.png')))))
 assert np.allclose(total.loc[d.Barcode],d.total_umi)
 data.append(dict(cancer='COAD',gene='UCKL1',sample=alias,patient=alias.split('_')[0],img=img,xy=d[['x','y']].to_numpy(),expr=v.loc[d.Barcode].to_numpy(),region=d.original_broad.to_numpy()))
def t(f,x,y,s,size=12,col=INK,bold=False,ha='left'):f.text(x,y,s,fontsize=size,color=col,weight='bold' if bold else 'normal',ha=ha,va='center')
def line(f,y):f.add_artist(Line2D([.05,.95],[y,y],transform=f.transFigure,color=LINE,lw=.8))
def base(title,sub):
 f=plt.figure(figsize=(16,9),facecolor='white');ax=f.add_axes([.044,.882,.168,.1]);ax.imshow(plt.imread(AS/'sysu_logo.png'));ax.axis('off');t(f,.95,.93,'CAMP  /  组会汇报',10,M,ha='right');line(f,.883);t(f,.05,.828,title,26,G,True);t(f,.05,.776,sub,11.5,M);return f
def save(f,name):
 for ext in ['png','pdf','svg']:f.savefig(P/(name+'.'+ext),dpi=200,facecolor='white')
 plt.close(f)
for c,g,num in [('BRCA','ASNS','03'),('COAD','UCKL1','10'),('PRAD','SLC6A6','07')]:
 es=[e for e in data if e['cancer']==c]
 sub={'BRCA':'CTA · 两位患者各一张代表切片 · 沿用 LYPLA1 的切片与病理标注','COAD':'Valdeolivas 2024 · 盲肠 / 右半结肠两病例 · 作者病理标注','PRAD':'Erickson 2022 · 同一患者两张切片 · 沿用 LYPLA1 的切片与病理标注'}[c]
 f=base(f'{c}：{g} 的空间表达与病理区域',sub)
 for xx,lab in zip([.16,.42,.68],['H&E','病理区域',g+' 实测表达']):t(f,xx,.720,lab,14,G,True)
 for e,y in zip(es,[.421,.094]):
  for j,xx in enumerate([.16,.42,.68]):
   ax=f.add_axes([xx,y,.225,.286]);img=e['img'];xy=e['xy'];v=e['expr'];reg=e['region'];ax.imshow(img,alpha=1 if j==0 else .52);ax.set_xlim(0,img.shape[1]);ax.set_ylim(img.shape[0],0);ax.axis('off')
   if j==1:
    for k in sorted(set(reg)):
     m=reg==k;ax.scatter(xy[m,0],xy[m,1],s=3.2,c=COL[k],lw=0,alpha=.88,rasterized=True)
   if j==2:
    ax.scatter(xy[:,0],xy[:,1],s=3.2,c='#E5E8EC',lw=0,rasterized=True);ix=np.argsort(v);ix=ix[v[ix]>0];ax.scatter(xy[ix,0],xy[ix,1],s=3.2,c=v[ix],cmap='viridis',vmin=0,vmax=CAP,lw=0,rasterized=True)
  t(f,.05,y+.259,e['sample'].replace('V19T26-','').replace('_Rep1',''),10.3,G,True);t(f,.05,y+.226,e['patient'],10,M)
  for k,rg in enumerate([z for z in COL if z in set(e['region'])]):
   yy=y+.179-k*.027;f.add_artist(Line2D([.057],[yy],transform=f.transFigure,marker='o',ls='',color=COL[rg],markersize=4.8));t(f,.068,yy,LABEL[rg],8.8)
  checks.append(dict(cancer=c,gene=g,section=e['sample'],spots=len(e['expr']),positive=int((e['expr']>0).sum()),clipped=int((e['expr']>CAP).sum())))
  for rg in sorted(set(e['region'])):
   vv=e['expr'][e['region']==rg];summ.append(dict(cancer=c,gene=g,section=e['sample'],patient=e['patient'],region=rg,spots=len(vv),mean=float(vv.mean()),detection=float((vv>0).mean())))
 # Separator stops before the shared colorbar.
 f.add_artist(Line2D([.05,.905],[.400,.400],transform=f.transFigure,color=LINE,lw=.8))
 cb=f.add_axes([.925,.29,.008,.27]);f.colorbar(ScalarMappable(norm=Normalize(0,CAP),cmap='viridis'),cax=cb,extend='max');cb.tick_params(labelsize=9);cb.set_ylabel(g+' · log1p(CP10K)',fontsize=9,color=M)
 line(f,.076);foot={'BRCA':'癌区与免疫区定位展示；未标注区不能当作正常上皮。','COAD':'S1 无合格非肿瘤上皮参照；S2 保留癌区与非肿瘤上皮，未新增 spot 级检验。','PRAD':'两张切片均来自同一患者；定位展示不等于多患者重复。'}[c]
 t(f,.05,.049,foot,10,G,True);t(f,.05,.024,'实测表达，无插补或平滑；统一色标 0–3，浅灰为零表达。',8.5,M)
 save(f,num+'_'+c+'_'+g+'_空间定位')
 print('SPATIAL',c,g,flush=True)
# GeoMx: retain ROI technology, never pretend this is a full tissue image.
d=B/'data/candidates/PDAC_GeoMx_Bell2025';m=pd.read_csv(d/'metadata_with_VI_subtypes_and_NGS.csv',index_col=0).set_index('Name');x=pd.read_csv(d/'ProbeQC_merged_batches.csv',index_col=0).loc[:,m.index];vst=pd.read_csv(d/'fully_batch_corrected_vsd.csv',index_col=0)
neg=x.loc['NegProbe-WTX'];loq=float(np.exp(np.log(neg).mean()+2*np.log(neg).std(ddof=1)));keep=(x.gt(loq).sum(1)>10)&(~x.index.str.startswith('NegProbe'));qc=dict(gene='SLC6A6',retained=bool(keep.get('SLC6A6',False)),loq=loq,above_loq=int((x.loc['SLC6A6']>loq).sum()))
(P/'PDAC_SLC6A6_GeoMx_QC.json').write_text(json.dumps(qc,indent=2))
f=base('PDAC：SLC6A6 的空间区域表达','Bell 2025 · GeoMx 区域测量 · 正常导管与癌上皮 ROI · 患者配对展示')
if qc['retained']:
 gx=x.loc[keep];q3=gx.quantile(.75,axis=0);ex=np.log2(x.loc['SLC6A6']/q3*np.exp(np.log(q3).mean())+1);roi=m[['Patient.Alias','Type']].copy();roi['q3_log2']=ex;roi['cpm_log2']=np.log2(x.loc['SLC6A6']/gx.sum(0)*1e6+1);roi['author_vst']=vst.loc['SLC6A6'] if 'SLC6A6' in vst.index else np.nan;roi.to_csv(R/'private_SLC6A6_geomx.tsv',sep='\t')
 summaries=[]
 for typ in ['ND','PDAC','VI','PNI']:
  ss=roi[roi.Type==typ];summaries.append(dict(region=typ,roi=len(ss),patients=ss['Patient.Alias'].nunique(),above_loq=int((x.loc['SLC6A6',ss.index]>loq).sum())))
 pd.DataFrame(summaries).to_csv(P/'PDAC_GeoMx_coverage.tsv',sep='\t',index=False)
 for k,(scale,label) in enumerate([('q3_log2','Q3 标准化'),('cpm_log2','CPM 敏感性'),('author_vst','作者 VST 敏感性')]):
  a=roi[roi.Type.isin(['ND','PDAC'])].groupby(['Patient.Alias','Type'])[scale].mean().unstack().dropna();a.to_csv(R/('private_geomx_pairs_'+scale+'.tsv'),sep='\t');delta=a.PDAC-a.ND;ax=f.add_axes([.10+k*.30,.29,.20,.37])
  for ref,test in a[['ND','PDAC']].to_numpy():ax.plot([0,1],[ref,test],c='#C77462' if test>ref else '#468B9E',lw=1.1,alpha=.65)
  for j,col in enumerate(['#468B9E','#C77462']):ax.scatter(np.full(len(a),j),a[['ND','PDAC']].iloc[:,j],s=26,c=col,zorder=3)
  ax.set_xticks([0,1],['正常导管','癌上皮']);ax.set_ylabel(scale,fontsize=10,color=M);ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=10);ax.set_title(label,fontsize=14,color=G,pad=18)
  t(f,.20+k*.30,.215,f'{int((delta>0).sum())} / {len(delta)} 患者癌区更高',13,G,True,ha='center')
  summ.append(dict(cancer='PDAC',gene='SLC6A6',section='GeoMx',patient='aggregate',region=scale,spots=len(a),mean=float(delta.mean()),detection=np.nan))
 t(f,.05,.113,'按患者汇总 ROI 后配对；本轮描述性展示，未新增显著性检验。',13,G,True)
 t(f,.05,.069,'GeoMx 为选定区域测量，不是完整切片热图；三种尺度不能直接比较数值。',10,M)
else:
 t(f,.50,.48,'SLC6A6 未通过既定 GeoMx 表达覆盖门槛',23,G,True,ha='center');t(f,.50,.36,'不据低于门槛的信号绘制癌区富集结论',16,M,ha='center')
save(f,'08_PDAC_SLC6A6_GeoMx')
pd.DataFrame(summ).to_csv(P/'spatial_region_summary.tsv',sep='\t',index=False)
(P/'spatial_display_checks.json').write_text(json.dumps(checks,indent=2));(P/'spatial_source_manifest.json').write_text(json.dumps(sources,indent=2))
print('SPATIAL_READY',flush=True)
