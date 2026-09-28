"""LYPLA1 in author-annotated Hirz Slide-seqV2; raw data stay on server165.
Usage: python prad_lypla1_spatial_v1.py RUN_ROOT
Requires R extraction script first. No imputation, smoothing, or reannotation.
"""
import sys, json, hashlib, platform
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.stats import mannwhitneyu
import statsmodels.api as sm
import statsmodels
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

root = Path(sys.argv[1]); pub = root/'public'; pub.mkdir(exist_ok=True)
d = pd.read_csv(root/'private/bead_values.tsv.gz', sep='\t')
assert len(d)==276884 and d.bead.is_unique
assert d[['xcoord','ycoord','total_umi','LYPLA1_log1p_cp10k']].notna().all().all()
assert (d.total_umi>=100).all()
assert np.allclose(d.LYPLA1_log1p_cp10k,np.log1p(1e4*d.LYPLA1_umi/d.total_umi))
# Metadata from Supplementary Data 2; grade mapping is not inferred from suffixes.
meta=pd.read_excel(root/'source/supplementary_data2.xlsx',sheet_name='Cases for Slide-seq',header=1)
assert {'Patient ID','Patient'} <= set(meta.columns)
md=dict(zip(meta['Patient ID'],meta['Patient']))
assert set(md)==set(d.section.unique())
d['tissue']=d.section.map(md)
d['detected']=d.LYPLA1_umi>0
sections=['Tumor01','Tumor02','Tumor07','Tumor08']
t=d[d.section.isin(sections)].copy()
assert set(t.tissue)=={'Tumor (HG)','Tumor (LG)'}

def summary(x):
    pos=x[x.detected]
    return dict(n_beads=len(x),n_positive=len(pos),detection_fraction=x.detected.mean(),
                mean_log1p_cp10k=x.LYPLA1_log1p_cp10k.mean(),
                mean_positive_log1p_cp10k=pos.LYPLA1_log1p_cp10k.mean() if len(pos) else np.nan,
                pseudobulk_cpm=1e6*x.LYPLA1_umi.sum()/x.total_umi.sum(),median_umi=x.total_umi.median())

group_rows=[]
for g,x in t.groupby('cell2'):
    row=dict(cancer='PRAD',gene='LYPLA1',cohort='Hirz2023_SlideSeqV2',group=g,**summary(x))
    ms=[summary(y) for _,y in x.groupby('section')]
    for k in ['detection_fraction','mean_log1p_cp10k','mean_positive_log1p_cp10k','pseudobulk_cpm']:
        row['section_equal_'+k]=np.mean([a[k] for a in ms])
    group_rows.append(row)
pd.DataFrame(group_rows).to_csv(pub/'tumor_celltype_summary.tsv',sep='\t',index=False)

contrasts=[]
for section,x in t.groupby('section'):
    for comp in ['Epithelial_Luminal','Nonmalignant_epithelial','Fibroblasts']:
        bmask=x.cell2.str.startswith('Epithelial_') if comp=='Nonmalignant_epithelial' else x.cell2.eq(comp)
        a=x[x.cell2=='Tumor']; b=x[bmask]; z=pd.concat([a,b]).copy()
        sa,sb=summary(a),summary(b)
        row=dict(section=section,tissue=md[section],comparison=comp,n_tumor=len(a),n_comparator=len(b))
        for k in sa: row['tumor_'+k]=sa[k]; row['comparator_'+k]=sb[k]
        row['mean_log_ratio']=sa['mean_log1p_cp10k']/sb['mean_log1p_cp10k']
        row['pseudobulk_ratio']=sa['pseudobulk_cpm']/sb['pseudobulk_cpm']
        row['p_bead_exploratory']=mannwhitneyu(a.LYPLA1_log1p_cp10k,b.LYPLA1_log1p_cp10k,alternative='two-sided').pvalue
        # Library-size offset removes sequencing-depth differences in count rate.
        # Fixed spatial grids supply robust sandwich SEs; these are sensitivity
        # checks, NOT patient-level replication or a full spatial covariance model.
        design=sm.add_constant((z.cell2=='Tumor').astype(float))
        for grid in [4,8]:
            ij=[]
            for coord in ['xcoord','ycoord']:
                vals=(z[coord]-x[coord].min())/(x[coord].max()-x[coord].min())
                ij.append(np.minimum((vals*grid).astype(int),grid-1))
            block=ij[0]*grid+ij[1]
            fit=sm.GLM(z.LYPLA1_umi,design,family=sm.families.Poisson(),offset=np.log(z.total_umi)).fit(
                cov_type='cluster',cov_kwds={'groups':block,'use_correction':True})
            row[f'rate_ratio_grid{grid}']=float(np.exp(fit.params.iloc[1]))
            row[f'p_cluster_grid{grid}']=float(fit.pvalues.iloc[1])
            ci=np.exp(fit.conf_int().iloc[1]);row[f'ci_low_grid{grid}']=ci.iloc[0];row[f'ci_high_grid{grid}']=ci.iloc[1]
            row[f'n_blocks_grid{grid}']=block.nunique()
        contrasts.append(row)
c=pd.DataFrame(contrasts)
c.to_csv(root/'private/section_contrasts.tsv',sep='\t',index=False)
public=[]
for comp,x in c.groupby('comparison'):
    public.append(dict(cancer='PRAD',gene='LYPLA1',comparison=comp,n_tumor_patients=2,n_tumor_sections=4,
        sections_higher_mean=int((x.mean_log_ratio>1).sum()),sections_higher_detection=int((x.tumor_detection_fraction>x.comparator_detection_fraction).sum()),
        sections_higher_positive_mean=int((x.tumor_mean_positive_log1p_cp10k>x.comparator_mean_positive_log1p_cp10k).sum()),
        sections_p_bead_lt_005=int((x.p_bead_exploratory<.05).sum()),
        min_p_bead=x.p_bead_exploratory.min(),max_p_bead=x.p_bead_exploratory.max(),
        min_mean_log_ratio=x.mean_log_ratio.min(),max_mean_log_ratio=x.mean_log_ratio.max(),
        min_pseudobulk_ratio=x.pseudobulk_ratio.min(),max_pseudobulk_ratio=x.pseudobulk_ratio.max(),
        sections_p_grid4_lt_005=int((x.p_cluster_grid4<.05).sum()),sections_p_grid8_lt_005=int((x.p_cluster_grid8<.05).sum()),
        min_p_grid4=x.p_cluster_grid4.min(),max_p_grid4=x.p_cluster_grid4.max(),
        min_p_grid8=x.p_cluster_grid8.min(),max_p_grid8=x.p_cluster_grid8.max(),
        patient_population_p='NOT_EVALUABLE_n2'))
pd.DataFrame(public).to_csv(pub/'contrast_summary.tsv',sep='\t',index=False)
# Grade-level summaries aggregate both sections per cancer donor; private only.
grade=[]
for (tissue,g),x in t.groupby(['tissue','cell2']):grade.append(dict(tissue=tissue,group=g,**summary(x)))
pd.DataFrame(grade).to_csv(root/'private/grade_group_summary.tsv',sep='\t',index=False)

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
colors={'Tumor':'#b2182b','Epithelial_Luminal':'#2166ac','Other epithelium':'#80cdc1','Other':'#dddddd'}
cmap=LinearSegmentedColormap.from_list('lypla1',['#eeeeee','#f6bcc2','#d94868','#70002e'])
fig,axs=plt.subplots(4,3,figsize=(12,14),constrained_layout=True)
rng=np.random.default_rng(20260928)
cap=float(np.quantile(t.loc[t.detected,'LYPLA1_log1p_cp10k'],.99))
for i,s in enumerate(sections):
    x=t[t.section==s].copy();x=x.iloc[rng.permutation(len(x))]
    x['display_group']=np.where(x.cell2=='Tumor','Tumor',np.where(x.cell2=='Epithelial_Luminal','Epithelial_Luminal',np.where(x.cell2.str.startswith('Epithelial_'),'Other epithelium','Other')))
    for g in colors:
        y=x[x.display_group==g];axs[i,0].scatter(y.xcoord,y.ycoord,s=1.6,c=colors[g],linewidths=0,rasterized=True,label=g)
    h=axs[i,1].scatter(x.xcoord,x.ycoord,c=x.LYPLA1_log1p_cp10k,s=1.6,cmap=cmap,vmin=0,vmax=cap,linewidths=0,rasterized=True)
    axs[i,1].set_title('LYPLA1 per bead\nRandom draw order')
    # Grid display summarizes observed counts; it is not imputation.
    h2=axs[i,2].hexbin(x.xcoord,x.ycoord,C=x.LYPLA1_log1p_cp10k,reduce_C_function=np.mean,gridsize=45,mincnt=5,cmap=cmap,vmin=0,vmax=.30,linewidths=0,rasterized=True)
    axs[i,2].set_title('Local mean, including zeros')
    axs[i,0].set_title(f'{s} | {md[s]}\nAuthor RCTD annotation')
    for ax in axs[i]:
        ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
        ax.set_xlim(x.xcoord.min()-20,x.xcoord.max()+20);ax.set_ylim(x.ycoord.min()-20,x.ycoord.max()+20)
fig.colorbar(h,ax=axs[:,1].tolist(),shrink=.45,label='log1p(CP10K)',extend='max')
fig.colorbar(h2,ax=axs[:,2].tolist(),shrink=.45,label='Mean log1p(CP10K)',extend='max')
axs[0,0].legend(markerscale=5,fontsize=7,loc='upper left')
fig.suptitle('PRAD LYPLA1 | 4 tumor sections, 2 patients\nSlide-seqV2 beads; author-inferred cell types, not pure single cells',fontsize=14)
for ext in ['png','pdf']:fig.savefig(pub/f'LYPLA1_spatial_maps.{ext}',dpi=220)
plt.close(fig)

fig,axs=plt.subplots(1,3,figsize=(11,3.8),constrained_layout=True)
cols=['tumor_mean_log1p_cp10k','tumor_detection_fraction','tumor_mean_positive_log1p_cp10k']
labels=['All beads: mean log1p(CP10K)','LYPLA1 detection (%)','Positive beads: mean log1p(CP10K)']
for ax,k,label in zip(axs,cols,labels):
    x=c[c.comparison=='Epithelial_Luminal'];other=k.replace('tumor_','comparator_');mul=100 if 'fraction' in k else 1
    for _,r in x.iterrows():
        co='#b2182b' if r.tissue=='Tumor (HG)' else '#e69f00'
        ax.plot([0,1],[r[other]*mul,r[k]*mul],'-o',color=co,alpha=.8,lw=1.5,ms=5)
    ax.set_xticks([0,1],['Author luminal','Author tumor']);ax.set_ylabel(label);ax.set_xlim(-.25,1.25)
fig.suptitle('Within-section comparison | 4 sections from 2 cancer patients\nRed: high grade; orange: low grade. Lines are sections, not independent patients.',fontsize=11)
for ext in ['png','pdf']:fig.savefig(pub/f'LYPLA1_spatial_comparison.{ext}',dpi=240)
plt.close(fig)

spec=dict(version='prad_lypla1_spatial_v1',gene='LYPLA1',seed=20260928,
    normalization='log1p(10000 * gene UMI / all matrix gene UMI); CPM ratio uses summed counts',
    primary_comparison='Author Tumor versus Epithelial_Luminal in each tumor section',
    secondary_comparisons=['all author Epithelial_*','Fibroblasts'],
    annotation='Original author RCTD labels; no LYPLA1-driven relabeling. Original reference feature inclusion not audited.',
    qc='Reuse author retained beads (singlet/doublet-certain); verify >=100 UMI; no new outcome-based filtering',
    statistical_unit='Bead for exploratory rank tests; spatial grid clusters for rate-ratio sensitivity. Cancer donors n=2; no population test.',
    tests='12 two-sided bead Mann-Whitney tests (4 sections x 3 comparisons); 24 Poisson offset sensitivity tests (4 and 8 grids). Nominal P; no FDR gate.',
    source_independence='Separate Hirz2023 study; no donor identity crosswalk to prior PRAD24 cohort, so no verified non-overlap claim.',
    plotting={'point_cap_positive_p99':cap,'local_mean_hex_grid':45,'mincnt':5,'local_mean_vmax':.30},
    software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,statsmodels=statsmodels.__version__,matplotlib=matplotlib.__version__))
(pub/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2)+'\n')
validation=dict(status='DONE',evidence_level='EXPLORATORY',n_beads=len(d),n_tumor_beads=len(t),n_sections=12,n_tumor_sections=4,n_cancer_patients=2,
    exact_annotation_and_coordinate_join=True,unique_beads=True,raw_nonnegative_integer_verified_by_R=True,
    normalization_recalculation_passed=True,min_umi=int(d.total_umi.min()),
    limitations=['RCTD labels include certain doublets; not pathology-only segmentation','Sparse gene detection; no imputation','Only 2 cancer patients','Grid robust P is a sensitivity analysis, not patient-level inference','LYPLA1 inclusion in original annotation model not audited'],
    original_labels_reused=True,prior_results_modified=False)
(pub/'validation.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2)+'\n')
urls={'hirz_annotations.rds':'https://drive.google.com/file/d/1lL721_aXt4-pQIfZcoiUxyvrgnaDtGdO/view',
      'hirz_coordinates.rds':'https://drive.google.com/file/d/10iNhgRAUOVy4RjfMkRnVihdpzJaEk8cC/view',
      'hirz_counts.rds':'https://drive.google.com/file/d/18dqg_SHMNKOtRtzSatnkcjuiACdqTC-L/view',
      'supplementary_data2.xlsx':'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-023-36325-2/MediaObjects/41467_2023_36325_MOESM5_ESM.xlsx',
      'hirz_article.html':'https://www.nature.com/articles/s41467-023-36325-2'}
manifest=[]
for name,url in urls.items():
    f=root/'source'/name
    if f.exists():
        h=hashlib.sha256()
        with f.open('rb') as stream:
            for chunk in iter(lambda:stream.read(1048576),b''):h.update(chunk)
        manifest.append(dict(file=name,url=url,sha256=h.hexdigest(),bytes=f.stat().st_size,residency='server165'))
pd.DataFrame(manifest).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
print(pd.DataFrame(public).to_string(index=False))
print(pd.DataFrame(group_rows).query("group in ['Tumor','Epithelial_Luminal','Fibroblasts']").to_string(index=False))
