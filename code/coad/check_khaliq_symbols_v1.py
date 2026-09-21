"""Check absent target/known alternate symbols without changing any mapping."""
import json
import sys
import zipfile
from pathlib import Path
from run_khaliq35_v1 import inner_csv

data,run=map(Path,sys.argv[1:3])
z=zipfile.ZipFile(data/'celltypes.joined.zip')
sets=[]
for name in z.namelist():
    if name.endswith('counts.csv.zip'):
        q,f=inner_csv(z,name); next(f)
        symbols=[line.split(',',1)[0].strip('"') for line in f]
        sets.append(set(symbols));f.close();q.close()
assert len(sets)==6
union=set.union(*sets)
candidates={'GSTT2':['GSTT2'],'SLC6A17':['SLC6A17','NTT4','MRT48'],'UPP2':['UPP2','UPASE2','UDRPASE2','UP2']}
result=dict(all_six_feature_sets_identical=all(x==sets[0] for x in sets),features_each=[len(x) for x in sets],checked_alias_hits={g:sorted(set(v)&union) for g,v in candidates.items()},paralog_GSTT2B_present='GSTT2B' in union,mapping_changed=False,sources=['https://www.ncbi.nlm.nih.gov/gene/388662','https://www.ncbi.nlm.nih.gov/protein/NP_001128570.1','https://www.ncbi.nlm.nih.gov/gene/2953'])
(run/'symbol_check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
