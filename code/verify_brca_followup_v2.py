"""Independent algebra/metadata checks; server only. No sample-level export."""
import argparse,csv,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
from prepare_brca_geo_metadata import convert

ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
s=a.root/'data/candidates/camp_primary_tissue_multicancer';o=a.out;agg=o/'aggregate'
raw=(o/'source/GSE37751_quick.txt').read_text()
assert convert(raw).decode().splitlines()==(o/'source/geo_metadata_lines.txt').read_text().splitlines()
rows=list(csv.reader((o/'source/geo_metadata_lines.txt').read_text().splitlines(),delimiter='\t'))
ids=rows[0][1:];md=pd.DataFrame(index=ids)
for row in rows[1:]:
    k=row[1].split(': ',1)[0];md[k]=[v.split(': ',1)[1] for v in row[1:]]
mp=pd.read_csv(s/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv',dtype=str);mp=mp[mp.Dataset.eq('BRCA1')]
gsm=mp.RNAID.str.replace('.CEL.gz','',regex=False)
assert md.loc[gsm,'tissue type'].map({'Tumor':'Tumor','Non-tumor':'Normal'}).tolist()==mp.TN.tolist()
tm=mp[mp.TN.eq('Tumor')];gsm=tm.RNAID.str.replace('.CEL.gz','',regex=False)
er=md.loc[gsm,'estrogen receptor status'].map({'Positive':1.,'Negative':0.}).to_numpy()
rna=pd.read_csv(s/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed'/tm.RNAFile.iloc[0],index_col=0)
met=pd.read_excel(s/'processed_metabolomics'/tm.MetabFile.iloc[0],sheet_name=tm.MetabFile_sheet.iloc[0],index_col=0);met.columns=met.columns.astype(str)
adj=pd.read_csv(agg/'er_adjusted_all_relations.tsv',sep='\t');maxerr=0.
for _,r in adj[adj.p_value.notna()].iterrows():
    x=met.loc[r.metabolite,tm.MetabID].to_numpy(float);y=rna.loc[r.gene,tm.RNAID].to_numpy(float)
    ok=np.isfinite(x)&np.isfinite(y)&np.isfinite(er);design=np.column_stack([np.ones(sum(ok)),er[ok]])
    rx,ry=rankdata(x[ok]),rankdata(y[ok]);ex=rx-design@np.linalg.lstsq(design,rx,rcond=None)[0];ey=ry-design@np.linalg.lstsq(design,ry,rcond=None)[0]
    maxerr=max(maxerr,abs(np.corrcoef(ex,ey)[0,1]-r.partial_rank_rho))
assert maxerr<1e-12
ok=adj.p_value.notna();q=multipletests(adj.loc[ok,'p_value'],method='fdr_bh')[1];qerr=max(abs(q-adj.loc[ok,'q_value']))
assert qerr<1e-12 and int(ok.sum())==171
loo=pd.read_csv(agg/'robustness_all_relations.tsv',sep='\t');focus=loo[loo.original_q.lt(.05)]
assert (focus.loo_same_sign==focus.loo_n).all()
summary={'GEO_conversion_exact':True,'all_108_tissue_labels_agree':True,'ER_partial_rank_independent_OLS_max_error':maxerr,'ER_BH_max_error':qerr,'ER_evaluable':int(ok.sum()),'original_significant_all_LOO_same_direction':True,'patient_independence':'NOT_VERIFIED'}
(agg/'independent_validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
