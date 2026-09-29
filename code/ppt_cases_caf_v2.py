"""Server-only CAF focus correction; no new inferential testing."""
from pathlib import Path
import sys,json,shutil,hashlib
import numpy as np,pandas as pd
R=Path(sys.argv[1]);R.mkdir(exist_ok=True);P=R/'public';P.mkdir(exist_ok=True)
OLD=R.parent/'20260929T_representative_cases_v1';AS=R.parent/'20260929T071500Z_four_sc_v3'
z=np.load(OLD/'PDAC_private.npz');m=json.loads((OLD/'PDAC_meta.json').read_text())
dp=R.parents[2]/'PDAC/B/20260928T031629Z_gse202051_lypla1_v1/private/cell_values.tsv.gz'
d=pd.read_csv(dp,sep='\t',low_memory=False);d=d[d.treatment=='Untreated'];assert len(d)==len(z['expr'])
assert set(z['types'])>=set(['Cancer-associated fibroblast','Malignant epithelial'])
# Same exact untreated cell order as the v1 extraction; no re-embedding.
x=pd.DataFrame(dict(patient=d.pid.to_numpy(),group=z['types'],expression=z['expr']))
a=x.groupby(['patient','group']).expression.agg(['mean','size']).reset_index();a=a[a['size']>=20]
n=a[a.group=='Malignant epithelial'].set_index('patient')['mean'];t=a[a.group=='Cancer-associated fibroblast'].set_index('patient')['mean'];ids=sorted(set(n.index)&set(t.index));pairs=np.c_[n.loc[ids],t.loc[ids]];delta=pairs[:,1]-pairs[:,0]
pd.DataFrame(dict(patient=ids,malignant=n.loc[ids],CAF=t.loc[ids])).to_csv(R/'private_CAF_pairs.tsv',sep='\t',index=False)
m.update(target='Cancer-associated fibroblast',locator='CAFs',ref_label='恶性上皮',test_label='CAF',n_pairs=len(pairs),higher=int((delta>0).sum()),lower=int((delta<0).sum()),equal=int((delta==0).sum()),mean_difference=float(delta.mean()),p=None,status='CAF_VS_MALIGNANT_DESCRIPTIVE',note='同患者 CAF 与恶性上皮各 ≥20 个细胞核；比较细胞类型，不是 CAF 肿瘤—正常变化。',conclusion=f'SLC6A6 表达以 CAF 更突出；{int((delta>0).sum())}/{len(pairs)} 位患者 CAF 高于恶性上皮。')
np.savez_compressed(R/'PDAC_private.npz',xy=z['xy'],types=z['types'],expr=z['expr'],pairs=pairs)
(R/'PDAC_meta.json').write_text(json.dumps(m,ensure_ascii=False,indent=2))
pd.DataFrame([m]).to_csv(P/'PDAC_CAF_patient_summary.tsv',sep='\t',index=False)
for f in ['msyh.ttc','msyhbd.ttc','sysu_logo.png']:shutil.copy2(AS/f,R/f)
code=(OLD/'ppt_four_sc_render_v3.py').read_text().replace("enumerate(['BRCA','COAD','PDAC','PRAD'],11)","enumerate(['PDAC'],5)").replace('LYPLA1','SLC6A6')
code=code.replace("f'{c}：从全细胞表达定位到同患者证据'","'PDAC：SLC6A6 在 CAF 中表达更突出'")
code=code.replace("if p is None:conclusion=f\"同患者描述：{meta['higher']}/{meta['n_pairs']} 位目标上皮表达更高；尚未作显著性检验。\"","if p is None:conclusion=meta['conclusion']")
code=code.replace("'患者目标组表达更高'","'患者 CAF 表达更高'").replace("['参照','目标']","['恶性上皮','CAF']")
code=code.replace("name=f'{page}_{c}_","name=f'{page:02d}_{c}_")
exec(compile(code,str(R/'caf_render.py'),'exec'),{'__name__':'__main__'})
files=[OLD/'PDAC_private.npz',OLD/'PDAC_meta.json',dp,OLD/'ppt_four_sc_render_v3.py']
pd.DataFrame([dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in files]).to_csv(P/'CAF_source_manifest.tsv',sep='\t',index=False)
print('CAF_RESULT',len(pairs),int((delta>0).sum()),float(delta.mean()),flush=True)
