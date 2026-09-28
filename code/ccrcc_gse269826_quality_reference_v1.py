"""Historical paired-cohort QC reference; raw data and identities stay on server."""
from pathlib import Path
import pandas as pd,numpy as np,h5py,json
from scipy import sparse
r=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716');s=r/'data/candidates/gse159115_kyn_cell_source_v0.1';out=r/'results/collaborative/ccRCC/B/20260928T033500Z_lypla1_gse269826_v1'
a=pd.read_csv(s/'GSE159115_ccRCC_anno.csv.gz');b=pd.read_csv(s/'GSE159115_normal_anno.csv.gz');t=a[a.anno.eq('Tumor')].copy();n=b[b.anno.isin(['PT-A','PT-B','PT-C'])].copy();t['condition']='Tumor';n['condition']='Normal';d=pd.concat([t,n]);c=d.groupby(['patient','condition']).size().unstack().fillna(0);paired=c.index[(c.Tumor>=20)&(c.Normal>=20)];d=d[d.patient.isin(paired)]
# Check author pct_MT units independently from a full raw H5 source.
p=sorted((s/'h5').glob('*.h5'))[0];sample='SI_'+p.name.split('_SI_',1)[1].split('_filtered')[0];obs=pd.concat([a,b]).set_index('cell')
with h5py.File(p) as f:
 z=f['GRCh38'];genes=[v.decode() for v in z['gene_names'][:]];bars=[sample+'_'+v.decode() for v in z['barcodes'][:]];x=sparse.csc_matrix((z['data'][:],z['indices'][:],z['indptr'][:]),shape=(len(genes),len(bars)));j=np.flatnonzero(pd.Index(bars).isin(obs.index));x=x[:,j];mi=[i for i,g in enumerate(genes) if g.startswith('MT-')];ratio=np.asarray(x[mi].sum(0)).ravel()/np.asarray(x.sum(0)).ravel();reported=obs.loc[np.array(bars)[j],'pct_MT'].to_numpy();err=float(np.max(np.abs(ratio-reported)));assert err<1e-12
rows=[]
for k,z in d.groupby('condition'):rows.append(dict(cohort='GSE159115',condition=k,n_patients=z.patient.nunique(),n_cells=len(z),median_UMI=z.total_UMI.median(),median_nFeature=z.no_genes.median(),median_mito_fraction=z.pct_MT.median(),mito_definition='author pct_MT is fraction, verified using raw counts'))
pd.DataFrame(rows).to_csv(out/'public/historical_cell_quality_summary.tsv',sep='\t',index=False)
(out/'public/historical_QC_unit_validation.json').write_text(json.dumps(dict(status='PASS',source_file=str(p),cells_checked=len(j),max_fraction_error=err,source_field='pct_MT',unit='fraction_not_percent'),indent=2))
