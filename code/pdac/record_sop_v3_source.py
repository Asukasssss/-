"""Record source-only scope after real numerical output and independent validation."""
import json
import pandas as pd
from record_sop_v3_resume import SOURCE,ROOT,VERSION,REPO,read,write,write_json,readme,stage,checksum

def main():
    s=json.loads((SOURCE/'summary.json').read_text());v=json.loads((SOURCE/'independent_validation.json').read_text());assert v['status']=='PASS' and s['genes']==687
    sp=json.loads((SOURCE/'analysis_spec.json').read_text());sp.update(status='DONE_SOURCE_NUMERICS',run_id=SOURCE.name,frozen_protocol_run_id='20260922T112300Z_sop_v3');write_json(SOURCE/'analysis_spec.json',sp)
    previous=read(ROOT/'06_EXTERNAL/20260922T112800Z_sop_v3_source/sc_dataset_registry.tsv')
    for cohort in s['cohorts']:
        idx=previous.cohort==cohort['cohort'];previous.loc[idx,'status']='DONE';previous.loc[idx,'expression_scale']='cellwise_log1p10k_then_equal_source_labels';previous.loc[idx,'reason']='Actual raw-count analysis complete;clinical cross-study identity not confirmed';previous.loc[idx,'n_cells']=cohort['cells'];previous.loc[idx,'n_source_labels']=cohort['source_labels'];previous.loc[idx,'n_genes_resolved']=cohort['matrix_resolved']
    write(SOURCE/'sc_dataset_registry.tsv',previous)
    cov=read(SOURCE/'sc_gene_coverage.tsv');top=read(SOURCE/'sc_source_stability.tsv')
    assert len(cov)==len(top)==2061 and not cov.duplicated(['cohort','stable_gene_id']).any()
    write_json(SOURCE/'validation.json',{'status':'PASS','historical_raw_matrix_hashes_match_before_run':True,'all687_in_all3cohorts':True,'full_gene_library_denominator':True,'cellwise_before_donor_aggregation':True,'donor_equal_means':True,'private_unit_tables_stayed_server':True,'source_no_P_q':True,'independent_verification':v,'figure_status':'PENDING_RENDER_QA','not_verified':['Protected clinical cross-study duplicate identities','CNV malignant confirmation','AmbientRNA correction','Disease differential expression','Metabolite relationship external validation','Mechanism']})
    stats='；'.join(f"{r['cohort']}：{r['cells']}细胞/{r['source_labels']}来源标签，{r['matrix_resolved']}/687基因可连接" for r in s['cohorts'])
    readme(SOURCE,'PDAC SOP v3：三队列全基因表达来源','全部687基因，包含250当前DIRECT及437条件/历史补充，不按患者P筛选。沿用作者肿瘤组织范围和预先固定粗细胞类别映射。',stats+'。每细胞按全基因库归一到1万后log1p，再按来源标签×细胞类型平均，最后等权汇总。覆盖阈值每标签类别20细胞、每类3标签；最低2可评估类别且最大检出比例≥1%才排名。1000次整标签重采样，80%标签仅描述排名稳定。','GSE242230的25个单位是作者样本标签，不能未经核实称25独立患者。三研究临床身份去重未完成。细胞最高来源不是唯一作用细胞，导管/上皮标签不确认恶性。公开图表仅显示群组汇总，私有逐供者/逐细胞数据留服务器。独立验证覆盖所有可评估汇总与排名，没有独立重读全部原始计数或重放全部bootstrap。','生成全687基因点图/热图并完成视觉核查，再接入关系/基因两表；不自动开展机制或细胞通讯。','服务器`python3 run_sop_v3_server.py --mode source --code-commit 321cd0f8f6ab8908a0d6daa190da7be57a8ed3a5`，随后`python3 verify_sop_v3_server.py --mode source`；图形`python code/pdac/plot_sop_v3_sources.py`。')
    stage(SOURCE,'06_EXTERNAL','687genes3cohorts cellwise sources;figures pending QA','Actual raw count/source stats independently checked;noSC P/q','Finish required dotplots and allgene heatmaps',status='PARTIAL')
    checksum(SOURCE);print(json.dumps(s))
if __name__=='__main__':main()
