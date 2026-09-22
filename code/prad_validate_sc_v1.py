"""Independent scalar source-profile and donor bootstrap checks, only aggregate output."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
def main(out):
 a=ad.read_h5ad(out/'source/PRAD24_cellxgene.h5ad',backed='r');obs=a.obs.copy();obs['ct']=obs.celltype_major_v2.astype(str);sel=obs.ct.eq('Epithelial');obs.loc[sel,'ct']='Epithelial_'+obs.loc[sel,'malignant_anno_merged'].astype(str)
 p=pd.read_csv(out/'private/sc_donor_profiles_private.tsv',sep='\t');groups=p[p.gene.eq('GSS')&p.n_cells.ge(20)].sample(6,random_state=4321);errors=[]
 col=int(np.flatnonzero(a.raw.var.feature_name.to_numpy()=='GSS')[0])
 for _,r in groups.iterrows():
  idx=np.flatnonzero(obs['type'].eq(r.partition)&obs.donor_id.eq(r.donor)&obs.ct.eq(r.celltype));x=a.raw.X[idx].tocsr();s=np.asarray(x.sum(1)).ravel();g=x[:,col].toarray().ravel();mu=np.log1p(g*10000/s).mean();det=np.mean(g>0);errors.append(abs(mu-r.mean_log1p_cp10k));assert abs(det-r.detection_fraction)<1e-12
 assert max(errors)<2e-5,max(errors)
 pp=p[p.partition.eq('cancer')];dn=sorted(pp.donor.unique());rng=np.random.default_rng(int(hashlib.sha256(('prad_sc_source_v1|cancer').encode()).hexdigest()[:8],16));draw=rng.integers(len(dn),size=(1000,len(dn)))
 z=pp[pp.gene.eq('GSS')&pp.n_cells.ge(20)&pp.celltype.ne('Unassigned')];eligible=z.groupby('celltype').donor.nunique();cts=eligible[eligible.ge(3)].index.tolist();means=z.groupby('celltype').mean_log1p_cp10k.mean().reindex(cts);top=means.idxmax();credits=[]
 for d in draw:
  labels=[dn[i] for i in d];scores={};detect={}
  for c in cts:
   v=z[z.celltype.eq(c)].set_index('donor');ls=[s for s in labels if s in v.index]
   if len(set(ls))>=3:scores[c]=float(v.loc[ls,'mean_log1p_cp10k'].mean());detect[c]=float(v.loc[ls,'detection_fraction'].mean())
  if len(scores)<2 or max(detect.values())<.01:continue
  best=max(scores.values());ties=[c for c,v in scores.items() if np.isclose(v,best,atol=1e-12,rtol=1e-12)];credits.append(1/len(ties) if top in ties else 0.)
 actual=pd.read_csv(out/'public/06_EXTERNAL/sc_source_stability.tsv',sep='\t');r=actual[actual.gene.eq('GSS')&actual.partition.eq('cancer')].iloc[0];assert r.top_celltype==top;assert abs(np.mean(credits)-r.bootstrap_top_frequency)<1e-12
 result=dict(status='DONE',scalar_mean_spotchecks=6,max_mean_error=float(max(errors)),detection_exact=True,independent_bootstrap_gene='GSS',independent_bootstrap_repeats=len(credits),bootstrap_frequency=float(np.mean(credits)),bootstrap_implementation_agrees=True)
 pub=out/'public/06_EXTERNAL';(pub/'independent_validation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);main(ap.parse_args().out)
