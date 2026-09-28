"""Package aggregate-only ccRCC LYPLA1 outputs and provenance on server165."""
import argparse,json,hashlib,shutil
from pathlib import Path
import pandas as pd

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for z in iter(lambda:f.read(1048576),b''):h.update(z)
 return h.hexdigest()
def main(out):
 r=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716');pub=out/'public';d=pd.read_csv(pub/'results.tsv',sep='\t');a=d[d['mode'].eq('primary')].iloc[0];qc=pd.read_csv(pub/'cell_quality_summary.tsv',sep='\t');oldpath=r/'results/collaborative/ccRCC/B/20260926T010000Z_lypla1_epithelial_v1/public/results.tsv';old=pd.read_csv(oldpath,sep='\t');old=old[old.comparator.eq('Normal_PT')].iloc[0]
 comp=pd.DataFrame([dict(cohort='GSE159115',n_pairs=old['n'],tumor_mean=old.tumor_mean,normal_mean=old.normal_mean,effect=old.effect,p_value=old.p_value,n_down=old.n_down,n_up=old.n_up),dict(cohort='GSE269826',n_pairs=a['n'],tumor_mean=a.tumor_mean,normal_mean=a.normal_mean,effect=a.effect,p_value=a.p_value,n_down=a.n_down,n_up=a.n_up)])
 comp.to_csv(pub/'cohort_comparison.tsv',sep='\t',index=False)
 paired=pd.read_csv(out/'private/paired_measurements.tsv',sep='\t');people=set(paired.loc[paired['mode'].eq('primary'),'patient']);cells=pd.read_csv(out/'private/cells_analysis.tsv',sep='\t');sel=json.loads((out/'source/ccrcc_gse269826_selection.json').read_text(encoding='utf-8-sig'));cells=cells[cells.patient.isin(people)];selected=cells[(cells.tissue.eq('ccRCC')&cells.annotation.isin(sel['tumor_labels']))|(cells.tissue.eq('Normal')&cells.annotation.isin(sel['normal_labels']))];qrows=[]
 for condition,z in selected.groupby('tissue'):
  qrows.append(dict(scope='primary_matched_patients',condition=condition,n_patients=z.patient.nunique(),n_cells=len(z),median_UMI=z.library.median(),median_nFeature=z.nFeature.median(),median_mito_fraction=z.mito.median(),mito_over_025_fraction=z.mito.gt(.25).mean(),doublet_flagged=int(z.doublet.ne('No').sum())))
 pd.DataFrame(qrows).to_csv(pub/'paired_cell_quality_summary.tsv',sep='\t',index=False)
 shutil.copyfile(out/'source/ccrcc_lypla1_cohort_selection.tsv',pub/'cohort_selection.tsv')
 dl=r/'data/candidates/ccrcc_lypla1_gse269826_20260928';shutil.copyfile(dl/'range_status.json',pub/'download_verification.json')
 (pub/'historical_status_correction.json').write_text(json.dumps(dict(previous='GSE269826 was described as incomplete based on residual aria2 markers',correction='Full Freshbiopsies object had been reconstructed and verified in previous Kynurenine work; current compressed SHA256 matches dbf7e0e5b9188aa1e1be8f8441a27ff0794248b2a511c51620bcca1d91b292f3. Refined epithelial object newly downloaded with validated ranges and gzip CRC. Old GSE159115 statistics unchanged.',historical_numbers_changed=False),indent=2))
 table='|分析|患者对数|癌减正常效应|P|q|下降/上升患者|\n|---|---:|---:|---:|---:|---:|\n'
 labels={'primary':'主分析','pseudobulk':'患者原始计数汇总','region_equal':'区域等权','min100cells':'每侧至少100细胞','strict_QC':'更严格质控','chr3p_loss':'限定作者3p缺失癌细胞'}
 for _,z in d.iterrows():table+=f"|{labels[z['mode']]}|{int(z['n'])}|{z.effect:.4f}|{z.p_value:.6g}|{z.q_value:.6g}|{int(z.n_down)}/{int(z.n_up)}|\n"
 direction='降低' if a.effect<0 else '升高';sig='达到预定主检验P<0.05' if a.p_value<.05 else '未达到预定主检验P<0.05'
 text=f'''# LYPLA1 更换单细胞队列后的癌上皮差异分析

## 本轮问题
使用更适合患者配对和癌上皮身份判定的新数据，比较 LYPLA1 在癌上皮与正常近端小管上皮的表达。数据选择和分析规则在查看目标基因结果前固定，不按显著性挑选队列。

## 输入与范围
选用 GSE269826（Lombardi等，Cell Reports 2025）新鲜活检对象。原研究新鲜对象132424细胞，原RNA计数36601基因。采用作者患者、病理、上皮细分及推断CNV注释；没有重聚类或自行猜测患者编号。排除培养、其他肾肿瘤病理和囊肿正常组织。最终{int(a['n'])}位患者配对，共癌上皮{int(a.tumor_cells)}、正常近端小管{int(a.normal_cells)}细胞；多区域按患者汇总。正常组织来自肿瘤邻近肾组织，并非健康志愿者。
该队列在可用配对患者数、恶性/近端小管注释以及QC可追溯性方面更适合本问题；不能据此声称其每项技术质控指标都优于旧数据。原作者允许较高线粒体比例以保留肾小管细胞，本轮另做更严格QC敏感性。

## 实际结果
主分析癌上皮均值{a.tumor_mean:.6f}，正常近端小管均值{a.normal_mean:.6f}，癌减正常{a.effect:.6f}；方向为{direction}，{sig}。检出率按配对患者等权：癌上皮{a.tumor_detection:.1%}、正常近端小管{a.normal_detection:.1%}。

{table}

主分析及除pseudobulk外的效应尺度是每细胞log1p(CP10K)的患者均值差，不能解释为倍数或直接转成表达下降百分比。pseudobulk为log2(1+CPM)差，不与其他行效应大小直接比较。主分析为用户指定单基因单比较，q=P；五项预定敏感性另行BH校正。每患者每条件至少20细胞。

## 新手解释
正常上皮参照固定为作者近端小管群，避免把全部肾上皮混合造成组成差异。统计单位是患者，不是细胞，也不是取样区域。旧GSE159115主比较为4对、P=0.375、癌减正常−0.091834；保留历史结果，不将两研究细胞合并计算显著性。当前研究不用于证明代谢酶活、通量或因果作用。

## 限制与反证
作者推断CNV支持癌细胞身份，但不是本轮新增基因型验证。邻近正常组织可能存在肿瘤相关背景变化。仅观察性独立队列复核，患者数仍有限。严格QC可能改变细胞组成，必须与主分析并列阅读。95%区间为10000次整患者配对bootstrap，和精确置换P不是同一检验的反演；以预定精确P判断显著性。留一患者结果见leave_one_patient_out_summary.tsv。
先前将GSE269826整体标记为未完成是依据残留下载标记的过度判断：此次已核对主对象完整校验及新下载的细分对象，见historical_status_correction.json。

## 当前决定
以新队列的预定主分析与全部敏感性共同解释，不选择性保留有利P值。原始数据、逐细胞和逐患者测量只保留在server165；公开文件仅为汇总及图。

## 下一步
若研究目标进一步涉及机制，需要蛋白/酶活或功能干预证据，不能由本轮RNA差异代替。

## 复现
先执行ccrcc_gse269826_extract_v1.R提取全基因分母和目标面板；完整下载后运行ccrcc_gse269826_refined_inspect_v1.R按条码对齐；按冻结selection执行export脚本；执行ccrcc_lypla1_gse269826_analyze_v1.py --out RUN；最后执行ccrcc_gse269826_validate_v1.R和finalize脚本。源码、参数及实际输入哈希见source_manifest.tsv。质量参考脚本仅复核旧队列QC，不重算旧效应。

来源：[GEO GSE269826](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269826)；[原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC7617885/)。
'''
 (pub/'README_CN.md').write_text(text,encoding='utf-8')
 # Manifest includes actual count-source caches and parent source checksums; no cell tables exported.
 manifest=[]
 for p in sorted((out/'source').glob('*'))+list((out/'private').glob('*cache.rds'))+[oldpath,r/'data/candidates/gse269826_ccrcc_v0.1/Freshbiopsies_inner.rds',dl/'epithelial_inner.rds']:
  if p.is_file():manifest.append(dict(source_id=p.name,path_or_url=str(p),sha256=sha(p),access_scope='SERVER_PRIVATE_DERIVED' if 'private' in p.parts or p.suffix=='.rds' else 'PUBLIC_CODE_OR_AGGREGATE'))
 download=json.loads((dl/'range_status.json').read_text());manifest.append(dict(source_id='GSE269826_Freshbiopsies_epithelial_reclustered.rds.gz',path_or_url=download['url'],sha256=download['sha256'],access_scope='SERVER_PRIVATE_MATRIX'))
 manifest.append(dict(source_id='GSE269826_Freshbiopsies.rds.gz',path_or_url=str(r/'data/candidates/gse269826_ccrcc_v0.1/GSE269826_Freshbiopsies.rds.gz'),sha256='dbf7e0e5b9188aa1e1be8f8441a27ff0794248b2a511c51620bcca1d91b292f3',access_scope='SERVER_PRIVATE_MATRIX_REHASHED_THIS_RUN'))
 for p in (r/'data/candidates/gse159115_kyn_cell_source_v0.1').glob('*anno.csv.gz'):manifest.append(dict(source_id=p.name,path_or_url=str(p),sha256=sha(p),access_scope='SERVER_PRIVATE_METADATA'))
 pd.DataFrame(manifest).to_csv(pub/'source_manifest.tsv',sep='\t',index=False)
 pd.DataFrame([dict(file=p.name,sha256=sha(p)) for p in sorted(pub.iterdir()) if p.is_file() and p.name!='checksums.tsv']).to_csv(pub/'checksums.tsv',sep='\t',index=False)
 (out/'.running').rename(out/'.done')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);main(p.parse_args().out)


