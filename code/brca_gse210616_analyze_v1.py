"""Exploratory marker association, not histology-defined malignant/normal validation."""
from pathlib import Path
import json,re,gzip,hashlib,sys,platform
import numpy as np,pandas as pd,h5py,scipy
from scipy import sparse,stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/data/candidates/BRCA_GSE210616')
OUT=Path(sys.argv[1]);OUT.mkdir(exist_ok=True,parents=True)
(OUT/'figures').mkdir(exist_ok=True)
PANELS={'Epithelial':['EPCAM','KRT8','KRT18','KRT19','KRT7','MUC1','KRT5','KRT14'],
        'Stromal':['COL1A1','COL1A2','COL3A1','DCN','LUM','COL6A1','PDGFRA','FAP'],
        'Immune':['PTPRC','LST1','TYROBP','CD3D','CD3E','CD79A','MS4A1'],
        'Endothelial':['PECAM1','VWF','KDR','EMCN']}
SPEC={'version':'v1','gene':'LYPLA1','status':'exploratory_marker_association','normalization':'log1p(count/total_counts*10000)',
      'QC':{'min_counts':500,'min_genes':200,'max_mito_fraction':0.25},'panels':PANELS,'score':'mean within-section gene z scores of log-normalized marker expression; LYPLA1 excluded',
      'region':'top panel score >0 and margin over second >=0.25; marker-inferred, not pathological annotation',
      'primary_requested_histology_test':'NOT_EVALUABLE_without_author_barcode_labels','unit':'patient; equal mean across their sections',
      'exploratory_tests':['depth-adjusted rank association LYPLA1 vs epithelial score','epithelial-marker-enriched vs stromal-marker-enriched mean log expression'],
      'min_region_spots':20,'sensitivity_min_spots':[10,30],'P':'two-sided exact sign test; nominal descriptive P; no spot-level P',
      'treatment':'NOT_EVALUABLE pending patient-specific clinical metadata','seed':20260927,
      'software':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'pandas':pd.__version__}}
(OUT/'analysis_spec.json').write_text(json.dumps(SPEC,indent=2))
def decode(a):return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in a])
def partial_rank(a,b,c):
    X=np.column_stack([np.ones(len(a)),stats.rankdata(c)])
    ra=stats.rankdata(a);rb=stats.rankdata(b)
    ra=ra-X@np.linalg.lstsq(X,ra,rcond=None)[0];rb=rb-X@np.linalg.lstsq(X,rb,rcond=None)[0]
    return np.corrcoef(ra,rb)[0,1]
rows=[];inputs=[];fails=[]
records=json.loads((ROOT/'manifest.json').read_text())
for rec in records:
    try:
        gsm=rec['gsm'];d=ROOT/gsm
        patient,section=map(int,re.search(r'patient (\d+), section (\d+)',rec['title']).groups())
        hp=next(d.glob('*.h5'))
        with h5py.File(hp) as f:
            g=f['matrix'];m=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=tuple(g['shape'][:]))
            genes=decode(g['features']['name'][:]);bars=decode(g['barcodes'][:])
        total=np.asarray(m.sum(axis=0)).ravel();ng=np.diff(m.indptr)
        mito=np.asarray(m[np.char.startswith(genes,'MT-')].sum(axis=0)).ravel()/np.maximum(total,1)
        ok=(total>=500)&(ng>=200)&(mito<=.25)
        m=m[:,ok];bars=bars[ok];total=total[ok]
        idx=np.flatnonzero(genes=='LYPLA1');assert len(idx)==1
        counts=m[idx[0]].toarray().ravel();y=np.log1p(counts/total*1e4)
        scores={};coverage={}
        for name,markers in PANELS.items():
            ii=np.flatnonzero(np.isin(genes,markers));coverage[name]=genes[ii].tolist();assert len(ii)>=3
            a=np.log1p(m[ii].toarray()/total[None,:]*1e4)
            sd=a.std(axis=1);a=(a-a.mean(axis=1)[:,None])/np.maximum(sd[:,None],1e-8)
            scores[name]=a.mean(axis=0)
        ss=np.column_stack(list(scores.values()));order=np.argsort(ss,axis=1)
        margin=ss[np.arange(len(y)),order[:,-1]]-ss[np.arange(len(y)),order[:,-2]]
        lab=np.array(list(scores),dtype=object)[order[:,-1]]
        lab[(ss.max(axis=1)<=0)|(margin<.25)]='Uncertain'
        ep=lab=='Epithelial';st=lab=='Stromal'
        row={'gsm':gsm,'patient':patient,'section':section,'n_input':len(ok),'n_qc':len(y),'n_detected':int((counts>0).sum()),
             'detected_fraction':float((counts>0).mean()),'n_epi':int(ep.sum()),'n_stroma':int(st.sum()),
             'rho_epi':stats.spearmanr(y,scores['Epithelial']).statistic,'partial_rho_epi':partial_rank(y,scores['Epithelial'],total),
             'rho_stroma':stats.spearmanr(y,scores['Stromal']).statistic,'partial_rho_stroma':partial_rank(y,scores['Stromal'],total)}
        for n in [10,20,30]:
            row[f'delta_{n}']=float(y[ep].mean()-y[st].mean()) if min(ep.sum(),st.sum())>=n else np.nan
        row['delta_pseudobulk']=float(np.log1p(counts[ep].sum()/total[ep].sum()*1e4)-np.log1p(counts[st].sum()/total[st].sum()*1e4)) if min(ep.sum(),st.sum())>=20 else np.nan
        pos=pd.read_csv(next(d.glob('*positions*.gz')),header=None,names=['barcode','in_tissue','array_row','array_col','pixel_row','pixel_col']).set_index('barcode')
        assert pos.index.is_unique and set(bars).issubset(pos.index)
        pos=pos.loc[bars];assert (pos.in_tissue==1).all()
        spot=pd.DataFrame({'barcode':bars,'LYPLA1_count':counts,'total_count':total,'LYPLA1_log1p':y,'marker_inferred_region':lab,**scores})
        spot.to_csv(OUT/f'{gsm}_spot_measurements.tsv.gz',sep='\t',index=False)
        scale=json.load(gzip.open(next(d.glob('*scalefactors*.gz')),'rt'))['tissue_hires_scalef']
        img=plt.imread(gzip.open(next(d.glob('*image*.gz')),'rb'),format='png')
        xx=pos.pixel_col.to_numpy()*scale;yy=pos.pixel_row.to_numpy()*scale
        fig,axs=plt.subplots(1,3,figsize=(15,5));colors={'Epithelial':'#e64980','Stromal':'#339af0','Immune':'#51cf66','Endothelial':'#f59f00','Uncertain':'#adb5bd'}
        for ax in axs:ax.imshow(img);ax.axis('off')
        axs[0].set_title('H&E')
        for label,color in colors.items():
            k=lab==label;axs[1].scatter(xx[k],yy[k],s=7,c=color,label=label,alpha=.85,linewidths=0)
        axs[1].set_title('Marker-inferred regions (not pathology)');axs[1].legend(loc='upper center',bbox_to_anchor=(.5,-.01),ncol=3,fontsize=7)
        vmax=max(float(np.quantile(y,.99)),.1)
        sc=axs[2].scatter(xx,yy,c=y,s=7,cmap='magma',vmin=0,vmax=vmax,linewidths=0)
        axs[2].set_title('LYPLA1 log-normalized expression');fig.colorbar(sc,ax=axs[2],shrink=.6)
        fig.suptitle(f'GSE210616 | Patient {patient}, section {section} | {gsm} | {len(y)} QC spots')
        fig.savefig(OUT/'figures'/f'{gsm}.png',dpi=150,bbox_inches='tight');plt.close(fig)
        rows.append(row)
        inputs.extend({'gsm':gsm,**f} for f in rec['files'])
        print(gsm,'ANALYZED',len(y),flush=True)
    except Exception as e:
        fails.append({'gsm':rec['gsm'],'error':repr(e)});print('ERROR',rec['gsm'],repr(e),flush=True)
sec=pd.DataFrame(rows);sec.to_csv(OUT/'section_results.tsv',sep='\t',index=False)
measures=['partial_rho_epi','rho_epi','partial_rho_stroma','delta_10','delta_20','delta_30','delta_pseudobulk']
pat=sec.groupby('patient')[measures].mean();pat.to_csv(OUT/'patient_results.tsv',sep='\t')
summary=[]
for col in measures:
    a=pat[col].dropna();nonzero=a[a!=0];n=len(nonzero);k=int((nonzero>0).sum())
    summary.append({'stage_id':'06_EXTERNAL','cancer':'BRCA','cohort':'GSE210616','gene':'LYPLA1','measure':col,'unit':'patient','n':len(a),'positive':k,'negative':int((a<0).sum()),'mean':float(a.mean()),'median':float(a.median()),'P_sign_two_sided':stats.binomtest(k,n,.5).pvalue if n else None,'status':'DONE','interpretation':'exploratory marker-based; not malignant versus normal'})
pd.DataFrame(summary).to_csv(OUT/'results.tsv',sep='\t',index=False)
pd.DataFrame(inputs).to_csv(OUT/'source_manifest.tsv',sep='\t',index=False)
validation={'sections_expected':43,'sections_analyzed':len(sec),'patients':len(pat),'spots_before_QC':int(sec.n_input.sum()),'spots_after_QC':int(sec.n_qc.sum()),'detected_spots':int(sec.n_detected.sum()),'failures':fails,'author_pathology_labels':'NOT_AVAILABLE_IN_DOWNLOADED_ESSENTIAL_FILES','treatment_stratification':'NOT_EVALUABLE','matrix_coordinate_joins':'validated per section','LYPLA1_in_annotation_markers':False,'patient_mapping':'explicit GEO sample titles','source_hashes':'SHA256 computed after successful downloads; no provider checksum comparison'}
(OUT/'validation.json').write_text(json.dumps(validation,indent=2))
fig,axs=plt.subplots(1,2,figsize=(10,5))
for ax,col,title in zip(axs,['partial_rho_epi','delta_20'],['Depth-adjusted association with epithelial markers','Epithelial vs stromal marker-enriched spots']):
    a=pat[col].dropna().sort_values();ax.barh([str(i) for i in a.index],a,color=['#d9485f' if v>0 else '#3b82b5' for v in a]);ax.axvline(0,color='black',lw=.8);ax.set_title(title,fontsize=10);ax.set_ylabel('Patient');ax.set_xlabel('Partial rank correlation' if col.startswith('partial') else 'Difference in mean log1p expression')
fig.tight_layout();fig.savefig(OUT/'figures'/'patient_summary.png',dpi=160);plt.close(fig)
print(json.dumps(validation),flush=True)
