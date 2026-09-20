"""Server-only PDAC source check and exact historical aggregate recovery."""
import argparse,csv,hashlib,json,platform,tarfile,traceback
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr

ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def read(p):
    with p.open(encoding='utf-8-sig',newline='') as f:
        d=csv.DictReader(f,delimiter='\t');return d.fieldnames,list(d)
def table(p,fields,rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)

def main():
    a=argparse.ArgumentParser();a.add_argument('--code-commit',required=True);args=a.parse_args()
    run=Path(__file__).resolve().parent
    assert run.parent==ROOT/'results/collaborative/PDAC/B'
    lock=run/'.running'
    with lock.open('x') as f:f.write('PDAC source_readiness.py\n')
    out=run/'public';out.mkdir(exist_ok=False)
    try:
        src=ROOT/'data/candidates/camp_primary_tissue_multicancer'
        mp=src/'metadata/MasterMapping_MetImmune_03_16_2022_release.csv'
        xp=src/'processed_metabolomics/PreprocessedData_PDAC.xlsx'
        rp=src/'gene_batch_reuse_20260910/pancancer_metabolomics/data/transcriptomics_processed/GSE62452.hugene10st.gene_symbol.csv'
        ap=ROOT/'results/CAMP_first_analysis_20260910/tables/CAMP_direct_metabolite_gene_associations.tsv'
        pp=ROOT/'results/CAMP_first_analysis_20260910/tables/CAMP_direct_gene_source_provenance.tsv'
        fp=run/'cancer_effects.tsv'
        assert sha(fp)=='1bb1d57480a0fbc2185d11f7598e67e7443aef9e40ea8a036c6da4e5503a4597'
        m=pd.read_csv(mp,dtype=str);m=m[m.Dataset=='PDAC'];tum=m[m.TN=='Tumor'];normal=m[m.TN=='Normal']
        rna=pd.read_csv(rp,index_col=0)
        met=pd.read_excel(xp,sheet_name='metabo_imputed_filtered_Tumor',index_col=0)
        norm=pd.read_excel(xp,sheet_name='metabo_imputed_filtered_Normal',index_col=0)
        raw=pd.read_excel(xp,sheet_name='data',index_col=0)
        _,frozen=read(fp);frozen=[r for r in frozen if r['cancer']=='PDAC'];sig=[r for r in frozen if float(r['effect_fdr'])<.05]
        coverage=[]
        for r in frozen:
            feature=r['feature_name'];count=list(met.index).count(feature)
            x=pd.to_numeric(met.loc[feature,tum.MetabID],errors='coerce') if count==1 else None
            rc=list(raw.index).count(feature)
            available=pd.to_numeric(raw.loc[feature,tum.MetabID],errors='coerce') if rc==1 else None
            coverage.append({'cancer':'PDAC','cohort':'PDAC','feature_name':feature,'metabolite_key':r['metabolite_key'],
              'original_q':r['effect_fdr'],'original_significant':float(r['effect_fdr'])<.05,
              'tumor_feature_matches':count,'normal_feature_matches':list(norm.index).count(feature),
              'n_tumor_processed_finite':int(np.isfinite(x).sum()) if x is not None else 'NA',
              'n_tumor_preimputation_available':int(np.isfinite(available).sum()) if available is not None else 'NA',
              'status':'DONE' if count==1 and rc==1 else 'NEEDS_REVIEW',
              'reason':'Source availability only; not verified mass-spectrometry detection'})
        table(out/'feature_coverage.tsv',list(coverage[0]),coverage)
        cols,old=read(ap);old=[r for r in old if r['cancer']=='PDAC']
        table(out/'historical_associations.tsv',cols,old)
        assert read(out/'historical_associations.tsv')[1]==old
        provenance=pd.read_csv(pp,sep='\t').set_index('path').sha256.to_dict()
        hashes={str(p):{'current':sha(p),'historical':provenance.get(str(p)),
                      'equal':provenance.get(str(p))==sha(p)} for p in [mp,xp,rp]}
        checks=[]
        for row in old:
            f=row['metabolite_name'];g=row['gene']
            assert list(met.index).count(f)==list(rna.index).count(g)==1
            x=pd.to_numeric(met.loc[f,tum.MetabID],errors='coerce').to_numpy(float)
            y=pd.to_numeric(rna.loc[g,tum.RNAID],errors='coerce').to_numpy(float)
            available=np.isfinite(pd.to_numeric(raw.loc[f,tum.MetabID],errors='coerce').to_numpy(float))
            keep=np.isfinite(x)&np.isfinite(y)
            if row['analysis_type']=='author_data_available_sensitivity':keep &= available
            rho=float(spearmanr(x[keep],y[keep]).statistic)
            checks.append({'gene':g,'metabolite_name':f,'analysis_type':row['analysis_type'],
              'actual_n':int(keep.sum()),'n_matches':int(row['n'])==int(keep.sum()),
              'rho_matches':bool(np.isclose(float(row['rho']),rho,rtol=0,atol=1e-12)),
              'availability_uses_identical_sample_mask':bool(available.all()),
              'review_note':'Historical sensitivity uses different random seed despite identical input; keep immutable, do not count as separate evidence'})
        table(out/'historical_statistics_check.tsv',list(checks[0]),checks)
        sources=[mp,xp,rp,ap,pp,fp,Path(__file__)]
        recover={
          'historical_cell_source_gene_summary.tsv':'results/PDAC_metabolite_gene_atlas_v0.1/tables/pdac_metabolite_gene_master.tsv',
          'historical_cell_source_relations.tsv':'results/PDAC_metabolite_gene_atlas_v0.1/tables/pdac_metabolite_linked_gene_results.tsv',
          'historical_cptac_effects.tsv':'results/PDAC_mechanism_exploration_v0.1/tables/cptac_candidate_rna_protein_effects.tsv',
          'historical_hypothesis_lock.yaml':'results/PDAC_mechanism_exploration_v0.1/parameters/hypothesis_lock.yaml',
          'historical_cell_source_parameters.yaml':'results/PDAC_metabolite_gene_atlas_v0.1/parameters/analysis_parameters.yaml'}
        inventory=[]
        for dest,rel in recover.items():
            p=ROOT/rel;sources.append(p);(out/dest).write_bytes(p.read_bytes())
            entry={'source_path':rel,'export_file':dest,'status':'NEEDS_REVIEW','reason':'Exact historical aggregate recovery; scope, independence and interpretation require review'}
            if p.suffix=='.tsv':
                c,rr=read(p);entry['rows']=len(rr);entry['genes']=len(set(r['gene'] for r in rr))
            else:entry['rows']='NA';entry['genes']='NA'
            inventory.append(entry)
        table(out/'historical_inventory.tsv',list(inventory[0]),inventory)
        summary={'cancer':'PDAC','n_mapped_tumor_specimens':len(tum),'n_mapped_normal_specimens':len(normal),
          'mapping_keys_complete':not m[['CommonID','MetabID','RNAID']].isna().any().any(),
          'mapping_keys_unique':all(m[k].is_unique for k in ['CommonID','MetabID','RNAID']),
          'tumor_rna_ids_found':int(tum.RNAID.isin(rna.columns).sum()),'tumor_metab_ids_found':int(tum.MetabID.isin(met.columns).sum()),
          'normal_rna_ids_found':int(normal.RNAID.isin(rna.columns).sum()),'normal_metab_ids_found':int(normal.MetabID.isin(norm.columns).sum()),
          'rna_gene_rows':len(rna),'rna_gene_labels_unique':rna.index.is_unique,
          'frozen_features':len(frozen),'significant_features':len(sig),
          'significant_features_present_unique_in_tumor':sum(r['original_significant'] and r['tumor_feature_matches']==1 for r in coverage),
          'significant_features_available_all27_before_imputation':sum(r['original_significant'] and r['n_tumor_preimputation_available']==len(tum) for r in coverage),
          'old_primary_relationships':sum(r['analysis_type']=='author_processed_primary' for r in old),
          'old_primary_q_lt_0_05':sum(r['analysis_type']=='author_processed_primary' and float(r['q_value'])<.05 for r in old),
          'patient_independence':'NOT_SEPARATELY_VERIFIED','tumor_normal_pairing':'NOT_INFERRED',
          'historical_input_hashes_match':all(v['equal'] for v in hashes.values())}
        dump(out/'summary.json',summary)
        dump(out/'input_hash_comparison.json',hashes)
        assert all(r['n_matches'] and r['rho_matches'] for r in checks)
        assert summary['mapping_keys_complete'] and summary['mapping_keys_unique']
        assert summary['tumor_rna_ids_found']==summary['tumor_metab_ids_found']==27
        assert summary['normal_rna_ids_found']==summary['normal_metab_ids_found']==12
        dump(out/'analysis_spec.json',{'analysis_version':'pdac_source_readiness_v1','run_id':run.name,
          'code_base_commit':args.code_commit,'script_sha256':sha(__file__),'unit':'author_mapped_specimens',
          'parameters':{'cancer':'PDAC','matrix_transforms':'none','new_imputation':False,'new_statistical_tests':False},
          'test_family':'Historical locked_direct_panel_v1 primary/sensitivity each 2 tests; unchanged',
          'seed':'Historical seed retained per row; no new permutation/bootstrap',
          'versions':{'python':platform.python_version(),'pandas':pd.__version__,'numpy':np.__version__,'scipy':scipy.__version__}})
        table(out/'source_manifest.tsv',['path','sha256','bytes'],[{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size} for p in sources])
        dump(out/'validation.json',{'status':'PASS','exact_historical_fields':True,'mapping_and_matrix_join_counts':True,
          'historical_n_rho_verified':True,'historical_source_hashes_match':summary['historical_input_hashes_match'],
          'not_verified':['patient independence','tumor-normal pairing','full42 direct gene mapping','external cohort independence','functional causality']})
        (out/'README_CN.md').write_text(f'''# PDAC 来源核对与历史结果接回

## 本轮问题
解除服务器访问阻塞，确认 PDAC 的矩阵覆盖，并接回适用的历史汇总。

## 输入与范围
服务器作者映射、PDAC 代谢矩阵、GSE62452 作者 RNA 矩阵及既有 PDAC 汇总。患者级数据和完整样本映射仅留服务器。历史统计列原样保留；feature_coverage.tsv 是来源覆盖计数，不是新统计检验。

## 实际结果
作者映射 27 个肿瘤标本、12 个正常标本，代谢/RNA 标识全部找到，三个映射键均完整唯一。303 条冻结特征保留；42 条原显著特征中 {summary['significant_features_present_unique_in_tumor']} 条在肿瘤矩阵唯一匹配，{summary['significant_features_available_all27_before_imputation']} 条在全部27标本的填补前作者数据可用。
旧关联只有鸟苷—PNP、鸟苷—SLC29A2 两条主分析，均 n=27，原新关联 q 均为0.4805。历史4行（含敏感性）逐字段保留，n/rho 已由当前输入核对；历史来源哈希一致：{summary['historical_input_hashes_match']}。
另接回单细胞14基因/20条关系汇总及 CPTAC 汇总，标记 NEEDS_REVIEW，不作为本轮独立验证完成。

## 新手解释
27 是作者映射的肿瘤标本数量；独立患者及肿瘤正常配对仍需来源证实。两条旧关联未通过 FDR，不代表其余40余条代谢特征均无候选。RNA/蛋白表达变化不等于酶活性或因果机制。

## 限制/反证
旧敏感性分析使用与主分析完全相同的有效样本，却换了随机种子，造成 P/q/区间略有差异。保留历史值并明确警示；后续同输入同方法应精确复用主统计，不能把此差异称为稳健性证据。
填补前 data 表的可用性不自动等于真正质谱检出。历史直接/近直接/上游候选需分别复核，单细胞表达不能证明功能干预；外部队列独立性尚未核实。

## 当前决定
本批来源核对及历史接回 DONE；03_PATIENT 总阶段 PARTIAL，02_MAPPING/06_EXTERNAL 为 NEEDS_REVIEW。解除 ACCESS_BLOCKED，不能宣称全候选关联完成。

## 下一步
以42条原显著特征建立完整直接人类酶/转运体关系清单；核实患者单位并锁定新检验族。适用旧统计复用，扩展后单独计算新q；同步审阅已接回的功能/外部证据。

## 复现命令
创建新的 `results/collaborative/PDAC/B/<RUN_ID>/`，复制本脚本和冻结 cancer_effects.tsv 后执行 `python3 source_readiness.py --code-commit <BASE_COMMIT>`。使用独占 .running 锁，已有输出拒绝覆盖。
''',encoding='utf-8')
        with tarfile.open(run/'public_delivery.tar','w') as tar:
            for p in sorted(out.iterdir()):tar.add(p,arcname=p.name)
        dump(run/'DONE.json',summary);lock.unlink();print(json.dumps(summary))
    except Exception:
        (run/'FAILED.txt').write_text(traceback.format_exc());raise

if __name__=='__main__':main()
