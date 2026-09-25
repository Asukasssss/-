from pathlib import Path
import pandas as pd,numpy as np,itertools,json
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716/results/collaborative/ccRCC/B/20260926T010000Z_lypla1_epithelial_v1');p=r/'public'
for ext in ['png','pdf']:(p/('LYPLA1_paired_epithelium.'+ext)).rename(r/'private'/('LYPLA1_paired_epithelium.'+ext))
d=pd.read_csv(p/'results.tsv',sep='\t');pro=pd.read_csv(r/'private/donor_profiles.tsv',sep='\t');checks=[]
for _,row in d.iterrows():
 z=pro[pro.cells.ge(20)].pivot(index='patient',columns='group',values='expression')[['Tumor',row.comparator]].dropna();delta=z.Tumor-z[row.comparator]
 checks.append(dict(comparator=row.comparator,n_match=len(z)==row['n'],effect_match=bool(np.isclose(delta.mean(),row.effect)),wilcoxon_exact_p=float(stats.wilcoxon(delta,method='exact').pvalue),primary_p=float(row.p_value)))
fig,ax=plt.subplots(figsize=(7,4));x=np.arange(2);ax.bar(x-.16,d.normal_mean,width=.32,label='Normal',color='#709fb0');ax.bar(x+.16,d.tumor_mean,width=.32,label='Tumor',color='#c66b65');ax.set_xticks(x);ax.set_xticklabels(['Proximal tubule','Defined renal tubular epithelium']);ax.set_ylabel('Patient-equal mean log1p(CP10K)');ax.set_title('LYPLA1 | 4 matched patients\nExact paired permutation P: 0.375 / 0.125');ax.legend();fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(p/('LYPLA1_aggregate_epithelium.'+ext),dpi=180)
(p/'independent_check.json').write_text(json.dumps(checks,indent=2))
print(pd.read_csv(p/'group_summaries.tsv',sep='\t').to_string(index=False));print(checks)
