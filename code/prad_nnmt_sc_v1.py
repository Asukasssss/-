"""Targeted NNMT supplement, reusing the frozen PRAD donor-equal engine."""
import argparse, json, shutil
from pathlib import Path
import anndata as ad
import pandas as pd
import prad_sc_source_v1 as engine

ap=argparse.ArgumentParser()
ap.add_argument('--out',type=Path,required=True)
ap.add_argument('--source-run',type=Path,required=True)
ap.add_argument('--code-commit',required=True)
args=ap.parse_args()
out=args.out
out.mkdir(exist_ok=False)
for folder in ['source','private','public']:(out/folder).mkdir()
(out/'.running').write_text('prad_discovery_v1')
for name in ['PRAD24_cellxgene.h5ad','PRAD24_download.json']:
 (out/'source'/name).symlink_to(args.source_run/'source'/name)
a=ad.read_h5ad(out/'source/PRAD24_cellxgene.h5ad',backed='r')
v=a.raw.var
hit=v[v.feature_name.astype(str).eq('NNMT')]
assert len(hit)==1, 'NNMT must match exactly once'
ensembl=str(hit.index[0]);a.file.close()
pd.DataFrame([{'gene':'NNMT','stable_gene_id':ensembl}]).to_csv(out/'source/direct_relations.tsv',sep='\t',index=False)
engine.main(out,args.code_commit)
pub=out/'public/06_EXTERNAL'
for stem in ['source_heatmap','source_dotplot']:
 for ext in ['png','pdf']:
  (pub/f'all101_{stem}.{ext}').unlink()
d=pd.read_csv(pub/'sc_celltype_profiles.tsv',sep='\t')
c=d[d.partition.eq('cancer') & d.status.eq('DONE')].sort_values('effect')
import matplotlib.pyplot as plt
fig,ax=plt.subplots(figsize=(9,5))
ax.barh(c.celltype,c.effect,color='#3287a8')
ax.set_xlabel('Donor-equal mean log1p(CP10K)')
ax.set_title('NNMT in PRAD cancer tissue | author cell annotations')
for i,(_,r) in enumerate(c.iterrows()):
 ax.text(r.effect+.015,i,f'{r.mean_detection_fraction:.1%} detected; n={int(r.n)} donors',va='center',fontsize=9)
ax.set_xlim(0,max(c.effect)*1.65);ax.spines[['top','right']].set_visible(False)
fig.tight_layout()
for ext in ['png','pdf']:fig.savefig(pub/f'NNMT_celltype_expression.{ext}',dpi=180)
plt.close(fig)
spec=json.loads((pub/'analysis_spec.json').read_text())
spec.update(target_selection='User-requested NNMT; separate targeted supplement, not added to original101',stable_gene_id_semantics='Exact raw.var Ensembl identifier',engine_version=engine.V)
engine.js(pub/'analysis_spec.json',spec)
# Independent mean/detection aggregation checks from retained donor-level data.
p=pd.read_csv(out/'private/sc_donor_profiles_private.tsv',sep='\t')
errors=[]
for _,r in d[d.status.eq('DONE')].iterrows():
 z=p[p.partition.eq(r.partition)&p.celltype.eq(r.celltype)&p.n_cells.ge(20)]
 errors.extend([abs(z.mean_log1p_cp10k.mean()-r.effect),abs(z.detection_fraction.mean()-r.mean_detection_fraction)])
assert max(errors)<1e-10
val=json.loads((pub/'validation.json').read_text())
val.update(independent_aggregation_max_error=max(errors),target='NNMT',target_ensembl=ensembl)
engine.js(pub/'validation.json',val)
manifest=pd.read_csv(pub/'source_manifest.tsv',sep='\t')
manifest=pd.concat([manifest,pd.DataFrame([dict(source_id=Path(__file__).name,path_or_url=str(Path(__file__)),sha256=engine.sha(__file__))])],ignore_index=True)
engine.save(manifest,pub/'source_manifest.tsv')
(pub/'README_CN.md').write_text('''# NNMT 单细胞定向补充

本轮问题：NNMT 在 PRAD 单细胞中的表达来源。

输入与范围：沿用前轮 CELLxGENE 固定版本和作者注释；肿瘤组织39755细胞、24供者，癌旁另列；NNMT 为用户指定补充，不修改原101基因结果。

实际结果：见 sc_celltype_profiles.tsv 和 sc_source_stability.tsv；主图仅肿瘤组织。每供者每类至少20细胞，每类至少3供者，供者等权；1000次整供者重抽样。raw.X全基因库归一化后 log1p(CP10K)。

新手解释：柱长为平均表达；标注百分比为先在各供者计算再等权平均的检测比例；n为合格供者数。排名保持率为描述性稳定性，不是P值。

限制/反证：单项研究；癌旁来自同研究且供者重叠，不是独立验证；不进行肿瘤与癌旁差异显著性检验；表达不等于酶活、代谢通量或1-MNA来源。保留前轮数据集GEO链接冲突警告，按CELLxGENE固定版本定位。

当前决定：仅作定向细胞表达背景，不声明机制成立。

下一步：独立队列与功能证据仍未完成。

复现命令：python prad_nnmt_sc_v1.py --out NEW_RUN --source-run ORIGINAL_RUN --code-commit COMMIT；同目录需 prad_sc_source_v1.py。
''',encoding='utf-8')
engine.save(pd.DataFrame([dict(file=f.name,sha256=engine.sha(f)) for f in sorted(pub.iterdir()) if f.name!='checksums.tsv']),pub/'checksums.tsv')
(out/'.running').rename(out/'COMPLETE')
print(c[['celltype','n','effect','mean_detection_fraction']].to_string(index=False))
print((pub/'sc_source_stability.tsv').read_text())
