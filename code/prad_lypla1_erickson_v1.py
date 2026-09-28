"""LYPLA1 in Erickson2022 author-consensus annotated Visium regions.
Usage: python prad_lypla1_erickson_v1.py RUN_ROOT BASE_COMMIT
Run only on server165; individual spot/section statistics remain private.
"""
import sys,json,hashlib,platform
from pathlib import Path
import numpy as np,pandas as pd,h5py,scipy,statsmodels
from scipy.sparse import csc_matrix
from scipy.stats import mannwhitneyu
import statsmodels.api as sm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

root=Path(sys.argv[1]);src=root/'source/Patient1';out=root/'public';out.mkdir(exist_ok=True)
private=root/'private';private.mkdir(exist_ok=True)
manifest=json.loads((root/'source/download_manifest.json').read_text())
for r in manifest:assert hashlib.sha256((root/r['file']).read_bytes()).hexdigest()==r['sha256']
assert len(manifest)==35
sections=sorted(x.name for x in src.iterdir() if x.is_dir());assert len(sections)==7
dfs=[];qc=[]
for s in sections:
    p=src/s
    with h5py.File(p/'filtered_feature_bc_matrix.h5') as h:
        m=h['matrix'];names=m['features/name'][:].astype(str);ids=m['barcodes'][:].astype(str)
        mat=csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=m['shape'][:])
        assert set(m['features/feature_type'][:].astype(str))=={'Gene Expression'}
    assert len(set(ids))==len(ids) and np.sum(names=='LYPLA1')==1
    assert np.all(mat.data>=0) and np.all(mat.data==np.round(mat.data))
    total=np.asarray(mat.sum(axis=0)).ravel();ng=np.asarray((mat>0).sum(axis=0)).ravel()
    gene=mat[np.flatnonzero(names=='LYPLA1')[0]].toarray().ravel()
    a=pd.read_csv(p/(s+'_Final_Consensus_Annotations.csv')).set_index('Barcode')
    pos=pd.read_csv(p/'tissue_positions_list.csv',header=None,names=['barcode','in_tissue','array_row','array_col','pixel_row','pixel_col']).set_index('barcode')
    assert a.index.is_unique and pos.index.is_unique and set(ids)<=set(pos.index)
    d=pd.DataFrame({'total_umi':total,'n_genes':ng,'LYPLA1_umi':gene},index=ids).join(a,validate='one_to_one').join(pos,validate='one_to_one')
    labels=set(a.Final_Annotations.dropna())
    assert all(not l.startswith('GG') or l in ['GG1','GG2','GG3','GG4','GG4 Cribriform','GG5'] for l in labels)
    keep=(d.total_umi>500)&(d.in_tissue==1)&d.Final_Annotations.notna()&d.Final_Annotations.ne('Exclude')
    z=d[keep].copy();z['section']=s;z['barcode']=z.index
    z['region']=np.where(z.Final_Annotations.str.startswith('GG'),'Cancer',z.Final_Annotations)
    z['expression']=np.log1p(1e4*z.LYPLA1_umi/z.total_umi);z['detected']=z.LYPLA1_umi>0
    qc.append(dict(section=s,n_matrix=len(d),n_annotation=len(a),n_annotation_matched=int(d.Final_Annotations.notna().sum()),
        n_after_qc=len(z),n_cancer=int((z.region=='Cancer').sum()),n_benign=int((z.region=='Benign').sum()),
        n_stroma=int((z.region=='Stroma').sum()),median_umi=z.total_umi.median(),median_genes=z.n_genes.median(),
        min_umi=z.total_umi.min(),max_umi=z.total_umi.max(),unmatched_annotation=len(set(a.index)-set(ids))))
    dfs.append(z)
d=pd.concat(dfs,ignore_index=True);d.to_csv(private/'spot_values.tsv.gz',sep='\t',index=False)
q=pd.DataFrame(qc);q.to_csv(private/'section_qc.tsv',sep='\t',index=False)

def metrics(x):
    y=x[x.detected]
    return dict(n=len(x),n_positive=len(y),detection=x.detected.mean(),mean=x.expression.mean(),
        positive_mean=y.expression.mean() if len(y) else np.nan,cpm=1e6*x.LYPLA1_umi.sum()/x.total_umi.sum(),
        median_umi=x.total_umi.median())
rows=[]
for s,x in d.groupby('section'):
    for ref in ['Benign','Stroma']:
        a=x[x.region=='Cancer'];b=x[x.region==ref];row=dict(section=s,reference=ref,n_cancer=len(a),n_reference=len(b))
        if min(len(a),len(b))<20:
            row.update(status='NOT_EVALUABLE',reason='Fewer than 20 QC spots in one comparison group');rows.append(row);continue
        ma,mb=metrics(a),metrics(b)
        row.update(status='DONE',reason='Within one patient; exploratory spot-level and spatial-grid sensitivity tests')
        row.update({'cancer_'+k:v for k,v in ma.items()});row.update({'reference_'+k:v for k,v in mb.items()})
        row['rate_ratio']=ma['cpm']/mb['cpm'];row['mean_difference']=ma['mean']-mb['mean']
        row['p_spot']=mannwhitneyu(a.expression,b.expression,alternative='two-sided').pvalue
        z=pd.concat([a,b]);design=sm.add_constant((z.region=='Cancer').astype(float))
        for grid in [4,8]:
            ij=[]
            for coord in ['array_row','array_col']:
                v=(z[coord]-x[coord].min())/(x[coord].max()-x[coord].min());ij.append(np.minimum((v*grid).astype(int),grid-1))
            block=ij[0]*grid+ij[1]
            row[f'blocks_grid{grid}']=block.nunique()
            row[f'cancer_blocks_grid{grid}']=block[z.region=='Cancer'].nunique()
            if block.nunique()<6 or block[z.region=='Cancer'].nunique()<2 or block[z.region==ref].nunique()<2:
                row[f'p_grid{grid}']=np.nan;continue
            fit=sm.GLM(z.LYPLA1_umi,design,family=sm.families.Poisson(),offset=np.log(z.total_umi)).fit(cov_type='cluster',cov_kwds={'groups':block,'use_correction':True})
            rr=float(np.exp(fit.params.iloc[1]));assert np.isclose(rr,row['rate_ratio'])
            row[f'p_grid{grid}']=float(fit.pvalues.iloc[1]);ci=np.exp(fit.conf_int().iloc[1]);row[f'low_grid{grid}']=ci.iloc[0];row[f'high_grid{grid}']=ci.iloc[1]
        rows.append(row)
c=pd.DataFrame(rows);c.to_csv(private/'section_contrasts.tsv',sep='\t',index=False)
summaries=[];standards=[]
for ref,z in c.groupby('reference'):
    x=z[z.status=='DONE']
    v=dict(cancer='PRAD',gene='LYPLA1',reference=ref,n_patients=1,n_sections_screened=7,n_sections_evaluable=len(x),
        n_sections_not_evaluable=int((z.status!='DONE').sum()),sections_higher_mean=int((x.mean_difference>0).sum()),
        sections_higher_detection=int((x.cancer_detection>x.reference_detection).sum()),
        sections_higher_positive_mean=int((x.cancer_positive_mean>x.reference_positive_mean).sum()),
        rate_ratio_min=x.rate_ratio.min(),rate_ratio_max=x.rate_ratio.max(),
        sections_p_spot_lt005=int((x.p_spot<.05).sum()),p_spot_min=x.p_spot.min(),p_spot_max=x.p_spot.max())
    for grid in [4,8]:
        v[f'n_evaluable_grid{grid}']=int(x[f'p_grid{grid}'].notna().sum())
        v[f'n_significant_grid{grid}']=int((x[f'p_grid{grid}']<.05).sum())
        v[f'p_grid{grid}_min']=x[f'p_grid{grid}'].min();v[f'p_grid{grid}_max']=x[f'p_grid{grid}'].max()
    summaries.append(v)
    standards.append(dict(cancer='PRAD',cohort='Erickson2022_Visium_Patient1',stage_id='06_EXTERNAL',run_id=root.name,
        analysis_version='prad_lypla1_erickson_v1',analysis_type='pathology_region_enrichment',metabolite_key='NA',metabolite_name='NA',gene='LYPLA1',
        unit='section_summary_within_one_patient',n=len(x),n_reference=1,effect_type='section_equal_mean_count_rate_ratio',effect=x.rate_ratio.mean(),
        ci_lower='NA',ci_upper='NA',p_value='NA',q_value='NA',test_family='Two comparisons per section; rank and two grid sensitivities',family_n_evaluable=int((c.status=='DONE').sum()),
        status='DONE',reason='Descriptive regional summary; no patient population P; nominal per-section P summarized separately',source_id='doi:10.17632/svw96g68dv.5',reference=ref))
pd.DataFrame(summaries).to_csv(out/'contrast_summary.tsv',sep='\t',index=False)
pd.DataFrame(standards).to_csv(out/'results.tsv',sep='\t',index=False)
# Public cohort-level summaries restricted to primary-comparison eligible sections.
eligible=c[(c.reference=='Benign')&(c.status=='DONE')].section.tolist()
main=d[d.section.isin(eligible)]
g=[dict(group=g,**metrics(x)) for g,x in main.groupby('region')]
pd.DataFrame(g).to_csv(out/'region_summary.tsv',sep='\t',index=False)

plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
cols={'Benign':'#2166ac','Cancer':'#b2182b','Stroma':'#bdbdbd'}
fig,axs=plt.subplots(1,3,figsize=(11,3.8),constrained_layout=True)
cc=c[(c.reference=='Benign')&(c.status=='DONE')].sort_values('section')
for ax,k,title in zip(axs,['mean','detection','positive_mean'],['All spots: mean log1p(CP10K)','Detection (%)','Positive spots: mean log1p(CP10K)']):
    for i,(_,r) in enumerate(cc.iterrows()):
        mul=100 if k=='detection' else 1
        ax.plot([0,1],[r['reference_'+k]*mul,r['cancer_'+k]*mul],'-o',color=plt.cm.tab10(i),label=r.section)
    ax.set_xticks([0,1],['Benign glands','Cancer']);ax.set_ylabel(title);ax.set_xlim(-.2,1.2);ax.spines[['top','right']].set_visible(False)
axs[0].legend(fontsize=7)
fig.suptitle(f'LYPLA1 | pathology consensus | {len(eligible)} comparable sections, ONE patient\nLines are tissue sections, not independent patients',fontsize=12)
for ext in ['png','pdf']:fig.savefig(out/f'LYPLA1_pathology_comparison.{ext}',dpi=240)
plt.close(fig)

cmap=LinearSegmentedColormap.from_list('lypla1',['#eeeeee','#f9c0c6','#d03b61','#720027'])
cap=float(d.loc[d.detected,'expression'].quantile(.99));rng=np.random.default_rng(20260928)
# Plot all sections, not a selection based on LYPLA1 outcome.
for page,sl in enumerate([sections[:4],sections[4:]],1):
    fig,axs=plt.subplots(len(sl),3,figsize=(11,3.5*len(sl)),constrained_layout=True,squeeze=False)
    for i,s in enumerate(sl):
        x=d[d.section==s].copy();x=x.iloc[rng.permutation(len(x))];p=src/s
        scales=json.loads((p/'scalefactors_json.json').read_text());scale=scales['tissue_hires_scalef'];im=plt.imread(p/(s+'_tissue_hires_image.png'))
        xx=x.pixel_col*scale;yy=x.pixel_row*scale
        for ax in axs[i]:ax.imshow(im);ax.set_xticks([]);ax.set_yticks([])
        axs[i,0].set_title(s+' | H&E')
        for group,color in cols.items():
            mask=x.region==group;axs[i,1].scatter(xx[mask],yy[mask],s=5,color=color,alpha=.8,linewidths=0,rasterized=True,label=group)
        axs[i,1].set_title('Author pathology consensus')
        h=axs[i,2].scatter(xx,yy,c=x.expression,s=6,cmap=cmap,vmin=0,vmax=cap,alpha=.85,linewidths=0,rasterized=True)
        axs[i,2].set_title('LYPLA1 | log1p(CP10K)')
        if i==0:axs[i,1].legend(fontsize=7,markerscale=2,loc='upper left')
    fig.colorbar(h,ax=axs[:,2].tolist(),shrink=.4,extend='max')
    fig.suptitle('Erickson et al. 2022 | Visium Patient 1\nMeasured spot expression, no imputation | H&E: Mendeley v5, CC BY-NC 3.0',fontsize=11)
    for ext in ['png','pdf']:fig.savefig(out/f'LYPLA1_pathology_maps_{page}.{ext}',dpi=210)
    plt.close(fig)

pd.DataFrame(manifest).to_csv(out/'source_manifest.tsv',sep='\t',index=False)
spec=dict(version='prad_lypla1_erickson_v1',base_git_commit=sys.argv[2],source_doi='10.17632/svw96g68dv.5',gene='LYPLA1',seed=20260928,
    selection='All seven Patient1/Visium_with_annotation sections, selected before expression analysis',
    annotation='Unmodified author Final_Consensus_Annotations; GG* cancer; Benign benign luminal glands; Stroma separate secondary reference',
    no_cnv_reference_selection=True,qc='In-tissue; total UMI >500; author Exclude and unannotated removed; >=20 spots per comparison group',
    normalization='log1p(10000*LYPLA1 UMI/all gene UMI); pseudobulk CPM from sums',
    tests='Two-sided Mann-Whitney per section; Poisson log(total UMI) offset with 4x4 and 8x8 spatial-grid cluster robust SE; nominal P only; no q gate',
    statistical_unit='Spot-level exploratory tests; one patient, no cohort-level inference',
    grid_eligibility='At least 6 occupied blocks overall and 2 per region',
    plots='All seven sections; common positive-expression p99 cap; random draw order; no interpolation',plot_cap=cap,
    software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,statsmodels=statsmodels.__version__,h5py=h5py.__version__,matplotlib=matplotlib.__version__),
    code_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in root.glob('prad_lypla1_erickson*.py')})
(out/'analysis_spec.json').write_text(json.dumps(spec,indent=2)+'\n')
validation=dict(status='DONE',n_patients=1,n_sections=7,n_primary_evaluable=len(eligible),n_spots_after_qc=len(d),n_primary_spots=len(main),
    raw_integer_counts=True,unique_barcodes_within_sections=True,unique_LYPLA1_feature=True,all_35_publisher_sha256_verified=True,
    normalization_finite=bool(np.isfinite(d.expression).all()),glm_rate_equals_pseudobulk_ratio=True,
    umi_median_range=[float(q.median_umi.min()),float(q.median_umi.max())],
    limitations=['One patient only','55 micron Visium spots contain multiple cells','Spatial-grid SE sensitivity is not full spatial modeling','Pathology region enrichment not pure cancer-cell expression'],
    QC_details_residency='server165/private/section_qc.tsv')
(out/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
print(pd.DataFrame(summaries).to_string(index=False));print(pd.DataFrame(g).query("group in ['Benign','Cancer','Stroma']").to_string(index=False));print(q.to_string(index=False))
