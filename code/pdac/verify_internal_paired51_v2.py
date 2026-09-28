"""Independent source-level n/rho, signed-rank DP and fixed-family BH checks."""
import argparse,json
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
from scipy.stats import rankdata,spearmanr
from internal_paired51_v2 import MP,XP,RP,AP,EXPECTED,validate_sources
def dp_p(d):
 d=d[d!=0]
 if not len(d):return 1.
 ranks=np.rint(rankdata(abs(d))*2).astype(int);counts=Counter({0:1})
 for r in ranks:
  new=Counter(counts)
  for k,v in counts.items():new[k+int(r)]+=v
  counts=new
 obs=int(ranks[d>0].sum());total=int(ranks.sum());dist=abs(2*obs-total)
 return sum(v for k,v in counts.items() if abs(2*k-total)>=dist)/2**len(d)
def bh(p):
 a=np.where(np.isfinite(p),p,1.);o=np.argsort(a);sortedq=np.minimum.accumulate((a[o]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1];q=np.empty(len(a));q[o]=np.minimum(sortedq,1.);q[~np.isfinite(p)]=np.nan;return q
def main():
 ap=argparse.ArgumentParser();ap.add_argument('run');args=ap.parse_args();run=Path(args.run);pub=run/'public';validate_sources(EXPECTED,[MP,XP,RP,AP])
 m=pd.read_csv(MP,dtype=str);m=m[m.Dataset=='PDAC'].copy();m['GSE ID']=m.RNAID.str.extract(r'^(GSM\d+)_',expand=False);author=pd.read_excel(AP,dtype=str);j=m.merge(author,on='GSE ID',validate='one_to_one');t=j[(j.TN=='Tumor')&j['Paired Sample ID'].notna()].sort_values('Paired Sample ID');groups=j.dropna(subset=['Paired Sample ID']).groupby('Paired Sample ID');pairs=[g for _,g in groups if set(g.TN)=={'Tumor','Normal'}];assert len(pairs)==11
 tp=[g.loc[g.TN=='Tumor','RNAID'].item() for g in pairs];npair=[g.loc[g.TN=='Normal','RNAID'].item() for g in pairs]
 rna=pd.read_csv(RP,index_col=0);met=pd.read_excel(XP,sheet_name='metabo_imputed_filtered_Tumor',index_col=0);raw=pd.read_excel(XP,sheet_name='data',index_col=0);a=pd.read_csv(pub/'associations.tsv',sep='\t');e=pd.read_csv(pub/'paired_RNA.tsv',sep='\t');checked=0
 for _,r in a.iterrows():
  if pd.isna(r.rna_label):assert r.status=='NOT_EVALUABLE';continue
  x=met.loc[r.metabolite_name,t.MetabID].to_numpy(float);y=rna.loc[r.rna_label,t.RNAID].to_numpy(float);mask=np.isfinite(x)&np.isfinite(y)
  if r.analysis_type=='availability':mask&=np.isfinite(raw.loc[r.metabolite_name,t.MetabID].to_numpy(float))
  assert int(mask.sum())==r.n
  if r.status=='DONE':assert abs(spearmanr(x[mask],y[mask]).statistic-r.effect)<1e-12;assert -1<=r.ci_lower<=r.ci_upper<=1;checked+=1
 for _,r in e.iterrows():
  if pd.isna(r.rna_label):assert r.status=='NOT_EVALUABLE';continue
  x=rna.loc[r.rna_label,tp].to_numpy(float);y=rna.loc[r.rna_label,npair].to_numpy(float);d=(x-y)[np.isfinite(x)&np.isfinite(y)];assert len(d)==r.n
  if r.status=='DONE':assert abs(dp_p(d)-r.p_value)<1e-12;assert abs(d.mean()-r.effect)<1e-12;assert (int((d>0).sum()),int((d<0).sum()),int((d==0).sum()))==(r.n_up,r.n_down,r.n_equal)
 for _,d in a.groupby('test_family'):assert np.allclose(bh(d.p_value.to_numpy(float)),d.q_value.to_numpy(float),equal_nan=True,atol=1e-12)
 assert np.allclose(bh(e.p_value.to_numpy(float)),e.q_value.to_numpy(float),equal_nan=True,atol=1e-12)
 result={'status':'PASS','source_level_association_n_rho_checks':checked,'RNA_signed_rank_independent_DP_checks':int((e.status=='DONE').sum()),'fixed_BH_families_checked':5,'permutation_P_independent_replay':'NOT_RUN;production method retained','patient_rows_exported':False}
 (pub/'independent_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
