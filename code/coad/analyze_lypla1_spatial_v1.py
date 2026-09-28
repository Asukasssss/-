"""Server-only LYPLA1 spatial descriptive analysis; source data never exported.

GSE283052: all 12 original frozen CRC sections, explicit barcode joins.
No clustering, region inference, smoothing, imputation, or significance tests.
"""
import argparse, gzip, hashlib, json, platform, re
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.io import mmread
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
from PIL import Image

ROOT = Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
DATA = ROOT/'data/candidates/coad_lypla1_spatial_20260928'
GENES = ['LYPLA1', 'EPCAM', 'COL1A1', 'PTPRC']
ap = argparse.ArgumentParser()
ap.add_argument('--out', type=Path, required=True)
ap.add_argument('--commit', required=True)
a = ap.parse_args()
assert a.out.parent == ROOT/'results/collaborative/COAD/B'
assert (a.out/'.running').is_dir()
PUB = a.out/'public'
PUB.mkdir(exist_ok=True)
FIG = PUB/'figures'
FIG.mkdir(exist_ok=True)
PRIVATE = a.out/'private'
PRIVATE.mkdir(exist_ok=True)
softpath = DATA/'GSE283052_family.soft.gz'
soft = gzip.open(softpath, 'rt').read()
samples = []
urls = {'GSE283052_family.soft.gz': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE283nnn/GSE283052/soft/GSE283052_family.soft.gz'}
for line in soft.splitlines():
    if line.startswith('!Series_supplementary_file') or line.startswith('!Sample_supplementary_file'):
        u = line.split(' = ', 1)[1].replace('ftp://', 'https://', 1)
        urls[u.rsplit('/', 1)[1]] = u
for block in soft.split('^SAMPLE = ')[1:]:
    gsm = block.splitlines()[0].strip()
    sample = re.search(r'!Sample_title = CRC tumor tissue (P\d+)', block).group(1)
    samples.append((gsm, sample))
assert len(samples) == 12 and len(set(s for _, s in samples)) == 12
spec = dict(analysis_version='COAD_LYPLA1_spatial_v1', code_lock_commit=a.commit,
    cohort='GSE283052', samples=samples, selection='All 12 source CRC sections selected before expression inspection',
    gene='LYPLA1', ensembl_id='ENSG00000120992', context_genes=GENES[1:],
    context_usage='Separate expression reference panels only; not cell labels or malignancy inference',
    unit='Visium capture spot; multi-cell; section summaries are descriptive',
    inclusion='Source in_tissue=1 AND total Gene Expression count>0; all others retained in audit',
    transform='log1p(10000 * gene UMI / total Gene Expression UMI per spot)',
    colors='gray-orange-red; one shared 0-to-observed-maximum scale per gene across all 12 sections; no clipping',
    primary_plot='LYPLA1 normalized expression on original tissue coordinates; no UMAP',
    sensitivity_plot='Raw LYPLA1 count and total-UMI maps',
    smoothing=False, imputation=False, tests='NONE', p_value=None, q_value=None,
    tumor_regions='NOT_EVALUABLE: no machine-readable pathologist spot labels acquired',
    patient_independence='12 original source sample codes; no sample-level linkage to CAMP independently verified',
    excluded_hd='GSE280315 LYPLA1 probe marked included=False; no zero-expression inference; HD acquisition stopped at eligibility gate',
    versions=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__))
(PUB/'analysis_spec.json').write_text(json.dumps(spec, indent=2)+'\n')
features = pd.read_csv(DATA/'GSE283052_features.tsv.gz', sep='\t', header=None, names=['id','gene','type'])
barcodes = pd.read_csv(DATA/'GSE283052_barcodes.tsv.gz', sep='\t', header=None)[0].astype(str)
assert barcodes.is_unique
gene_rows = {}
for gene in GENES:
    ix = np.flatnonzero((features.gene == gene) & (features.type == 'Gene Expression'))
    assert len(ix) == 1, (gene, ix)
    gene_rows[gene] = int(ix[0])
assert features.iloc[gene_rows['LYPLA1']].id == 'ENSG00000120992'
gene_mask = (features.type == 'Gene Expression').to_numpy()
manifest_paths = [softpath, DATA/'GSE283052_features.tsv.gz', DATA/'GSE283052_barcodes.tsv.gz']
records, summary, checks = [], [], []
for gsm, sample in samples:
    prefix = gsm+'_'+sample
    files = {k:DATA/(prefix+suf) for k,suf in {
        'matrix':'matrix.mtx.gz', 'positions':'tissue_positions_list.csv.gz',
        'scale':'scalefactors_json.json.gz','image':'tissue_lowres_image.png.gz'}.items()}
    for p in files.values():
        assert p.is_file(), str(p)
    manifest_paths.extend(files.values())
    with gzip.open(files['matrix'], 'rb') as f, threadpool_limits(limits=1):
        matrix = mmread(f).tocsr()
    assert matrix.shape == (len(features), len(barcodes)), (sample, matrix.shape)
    assert np.isfinite(matrix.data).all() and (matrix.data >= 0).all()
    assert np.equal(matrix.data, np.floor(matrix.data)).all()
    total = np.asarray(matrix[gene_mask].sum(axis=0)).ravel().astype(float)
    n_genes = np.asarray((matrix[gene_mask] > 0).sum(axis=0)).ravel()
    raw = {g:matrix[i].toarray().ravel().astype(float) for g,i in gene_rows.items()}
    del matrix
    pos = pd.read_csv(files['positions'], header=None,
        names=['barcode','in_tissue','array_row','array_col','pixel_row','pixel_col'])
    assert pos.barcode.is_unique and pos.in_tissue.isin([0,1]).all()
    joined = pd.DataFrame({'barcode':barcodes, 'total_umi':total, 'n_genes':n_genes}).merge(
        pos, on='barcode', how='left', validate='one_to_one', indicator=True)
    assert (joined._merge == 'both').all()
    joined = joined.drop(columns='_merge')
    assert np.array_equal(joined.barcode.to_numpy(), barcodes.to_numpy())
    scale = json.load(gzip.open(files['scale'], 'rt'))
    with gzip.open(files['image'], 'rb') as f:
        im = np.array(Image.open(f).convert('RGB'))
    sf = float(scale['tissue_lowres_scalef'])
    joined['x'] = joined.pixel_col * sf
    joined['y'] = joined.pixel_row * sf
    eligible = (joined.in_tissue == 1) & (joined.total_umi > 0)
    assert eligible.any()
    assert np.isfinite(joined.loc[eligible, ['x','y']]).all().all()
    assert joined.loc[eligible,'x'].between(0, im.shape[1]).all()
    assert joined.loc[eligible,'y'].between(0, im.shape[0]).all()
    for gene in GENES:
        joined[gene+'_raw'] = raw[gene]
        joined[gene+'_log1pCP10K'] = np.where(total > 0, np.log1p(10000*raw[gene]/np.maximum(total,1)), np.nan)
    joined['eligible'] = eligible
    # All spot-level values remain exclusively on the server.
    joined.to_csv(PRIVATE/(sample+'_spot_values.tsv.gz'), sep='\t', index=False)
    r = joined.loc[eligible].copy()
    records.append(dict(sample=sample, gsm=gsm, frame=r, image=im, scale=scale,
        diameter=float(scale['spot_diameter_fullres'])*sf))
    for gene in GENES:
        v=r[gene+'_log1pCP10K']; c=r[gene+'_raw']
        summary.append(dict(cancer='COAD',cohort='GSE283052',stage_id='06_EXTERNAL',run_id=a.out.name,
            analysis_version=spec['analysis_version'],analysis_type='spatial_expression_description',
            metabolite_key='NA',metabolite_name='NA',gene=gene,unit='in_tissue_spot_with_positive_total_UMI',
            n=len(r),n_reference='NA',effect_type='mean_log1pCP10K',effect=float(v.mean()),
            ci_lower='NA',ci_upper='NA',p_value='NA',q_value='NA',test_family='NONE',family_n_evaluable=0,
            status='DONE',reason='descriptive only; spot mixtures; no author region labels acquired',source_id=gsm,
            sample=sample,ensembl_id=features.iloc[gene_rows[gene]].id,
            positive_spots=int((c>0).sum()),detected_fraction=float((c>0).mean()),
            median_log1pCP10K=float(v.median()),max_log1pCP10K=float(v.max()),
            total_gene_UMI=int(c.sum()),max_gene_UMI=int(c.max()),
            total_UMI=int(r.total_umi.sum()),median_total_UMI=float(r.total_umi.median()),
            median_detected_genes=float(r.n_genes.median()),pooled_CP10K=float(10000*c.sum()/r.total_umi.sum())))
    checks.append(dict(sample=sample,source_id=gsm,matrix_genes=len(features),matrix_barcodes=len(barcodes),
        all_barcodes_unique=True,all_count_spots_position_matched=True,coordinate_bounds_pass=True,
        source_in_tissue=int((joined.in_tissue==1).sum()),eligible_spots=int(eligible.sum()),
        in_tissue_zero_total=int(((joined.in_tissue==1)&(total==0)).sum()),
        out_of_tissue=int((joined.in_tissue==0).sum()),spot_diameter_lowres_pixels=float(scale['spot_diameter_fullres'])*sf,
        image_height=int(im.shape[0]),image_width=int(im.shape[1]),gene_count_sum=int(raw['LYPLA1'].sum())))
    print('ANALYZED',sample,len(r),int((r.LYPLA1_raw>0).sum()),flush=True)

pd.DataFrame(summary).to_csv(PUB/'expression_summary.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(PUB/'sample_qc.tsv',sep='\t',index=False)
cmap=LinearSegmentedColormap.from_list('gray_orange_red',['#dedede','#fdbb84','#ef6548','#b30000','#67000d'])
limits={g:max(float(r['frame'][g+'_log1pCP10K'].max()) for r in records) for g in GENES}
rawmax=max(float(r['frame'].LYPLA1_raw.max()) for r in records)
depthmax=max(float(np.log1p(r['frame'].total_umi).max()) for r in records)

def draw(ax,record,values,vmax,title,background=True):
    df=record['frame'];im=record['image']
    if background:ax.imshow(im,alpha=.70)
    circles=[Circle((x,y),radius=record['diameter']/2) for x,y in zip(df.x,df.y)]
    patches=PatchCollection(circles,cmap=cmap,norm=Normalize(0,vmax),linewidth=0,rasterized=True)
    patches.set_array(np.asarray(values,dtype=float));ax.add_collection(patches)
    ax.set_xlim(0,im.shape[1]);ax.set_ylim(im.shape[0],0);ax.set_aspect('equal');ax.axis('off');ax.set_title(title,fontsize=11)
    return patches

plt.rcParams.update({'font.size':10,'savefig.facecolor':'white'})
fig,axes=plt.subplots(3,4,figsize=(14,12))
for ax,r in zip(axes.flat,records):
    f=r['frame'];fraction=(f.LYPLA1_raw>0).mean()
    p=draw(ax,r,f.LYPLA1_log1pCP10K,limits['LYPLA1'],f"{r['sample']} | {len(f):,} spots | {fraction:.1%} detected")
fig.subplots_adjust(left=.025,right=.91,bottom=.08,top=.92,wspace=.08,hspace=.18)
fig.colorbar(p,cax=fig.add_axes([.935,.22,.015,.55]),label='LYPLA1 log1p(UMI / total UMI x 10,000)')
fig.suptitle('LYPLA1 | original CRC tissue locations | GSE283052',fontsize=17,y=.975)
fig.text(.03,.035,'All 12 sections; shared full-range scale; gray = zero detected count. Spots contain multiple cells.\nH&E background with expression overlay; no smoothing, region labeling or significance test.',fontsize=10)
fig.savefig(FIG/'LYPLA1_spatial_all12.png',dpi=160)
fig.savefig(FIG/'LYPLA1_spatial_all12.pdf')
plt.close(fig)
for r in records:
    f=r['frame'];fig,axes=plt.subplots(2,3,figsize=(14,10))
    panels=[('LYPLA1_log1pCP10K',limits['LYPLA1'],'LYPLA1 | log1p CP10K'),
        ('LYPLA1_raw',rawmax,'LYPLA1 | raw UMI'),('depth',depthmax,'Library size | log1p total UMI')]
    panels += [(g+'_log1pCP10K',limits[g],g+' | log1p CP10K') for g in GENES[1:]]
    for ax,(col,lim,title) in zip(axes.flat,panels):
        val=np.log1p(f.total_umi) if col=='depth' else f[col]
        p=draw(ax,r,val,lim,title)
        fig.colorbar(p,ax=ax,fraction=.034,pad=.02)
    fig.suptitle(r['sample']+' | LYPLA1 and expression context | GSE283052',fontsize=16)
    fig.text(.02,.018,'Gene-specific shared scales across 12 sections. EPCAM/COL1A1/PTPRC are expression references, not inferred cell types.\nMulti-cell capture spots; no malignant/non-malignant labels assigned; no tests or smoothing.',fontsize=10)
    fig.tight_layout(rect=[0,.065,1,.94])
    fig.savefig(FIG/(r['sample']+'_LYPLA1_context.png'),dpi=135)
    fig.savefig(FIG/(r['sample']+'_LYPLA1_context.pdf'))
    plt.close(fig)

# Audit excluded HD sources without pretending missing measurement is zero.
hd=[]
for gsm,sample in [('GSM8594567','P1CRC'),('GSM8594568','P2CRC'),('GSM8594569','P5CRC')]:
    p=DATA/(gsm+'_'+sample+'_probe_set.csv.gz')
    row=dict(cohort='GSE280315',source_id=gsm,sample=sample,gene='LYPLA1',status='NOT_EVALUABLE',
        reason='probe not eligible; full matrix acquisition not completed',probe_checked=False,included='NA')
    if p.exists():
        probes=pd.read_csv(p,comment='#')
        sub=probes.loc[probes.gene_id=='ENSG00000120992']
        assert len(sub)==1 and not bool(sub.iloc[0]['included'])
        row.update(probe_checked=True,included=False,probe_id=sub.iloc[0].probe_id)
        manifest_paths.append(p)
    else:row['reason']='sample probe file not acquired; cohort held after excluded probe observed in acquired samples'
    hd.append(row)
pd.DataFrame(hd).to_csv(PUB/'hd_eligibility.tsv',sep='\t',index=False)
manifest=[]
for path in manifest_paths:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    url=urls.get(path.name)
    if url is None:
        gsm=path.name.split('_')[0]
        url='https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8594nnn/'+gsm+'/suppl/'+path.name
    manifest.append(dict(source_id=path.name,url=url,server_path=str(path),bytes=path.stat().st_size,sha256=h.hexdigest()))
pd.DataFrame(manifest).to_csv(PUB/'source_manifest.tsv',sep='\t',index=False)
validation=dict(status='PASS',scope='12-section identity/count/coordinate/transform/figure-generation checks',
    samples=len(records),gene_identity_unique=True,all_matrix_dimensions_match=True,
    explicit_barcode_joins=True,all_eligible_coordinates_in_image=True,integer_nonnegative_counts=True,
    all12_same_LYPLA1_scale=True,normalized_scale_maxima=limits,raw_LYPLA1_max=rawmax,
    sum_spots=sum(c['eligible_spots'] for c in checks),p_values_computed=False,
    unverified=['Machine-readable histology region labels unavailable','Cross-study patient independence not audited',
        'No pathologist review of generated maps','Raw reads and original image alignment not independently reprocessed'])
(PUB/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
(a.out/'ANALYSIS_DONE').write_text('DONE\n')
print('ANALYSIS_DONE',validation['sum_spots'],flush=True)
