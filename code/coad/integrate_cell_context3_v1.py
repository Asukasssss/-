"""Append third-study context without altering any historical relation field."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

RUN='20260921T095621Z_cell_context3_v1'
OLD='results/COAD/06_EXTERNAL/20260920T091100Z_cell_paired35_v1/relations39_integrated.tsv'
LEE='results/COAD/06_EXTERNAL/20260920T081200Z_cell_expression35_v2/tumor_lineage_overview.tsv'
KBASE='results/COAD/06_EXTERNAL/20260921T085412Z_khaliq35_v1'
IBASE='results/COAD/06_EXTERNAL/20260921T082557Z_integral_admission_v1'
CORRECTION='https://link.springer.com/article/10.1186/s13059-022-02724-9'

def read(path):
    with path.open(encoding='utf-8',newline='') as f:
        r=csv.DictReader(f,delimiter='\t');return r.fieldnames,list(r)

def write(path,fields,rows):
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=a.repo;d=a.out;d.mkdir(parents=True,exist_ok=True)
    oldfields,old=read(r/OLD);_,historical=read(r/LEE);_,khaliq=read(r/KBASE/'tumor_lineage_overview.tsv');_,primary=read(r/KBASE/'tumor_lineage_all35.tsv')
    assert len(old)==39 and len({x['gene'] for x in old})==35
    assert hashlib.sha256((r/OLD).read_bytes()).hexdigest()=='db476e7d53ed2c1b86df01dd3c5562830e5836424077b0e788035be6fa89aba7'
    assert len({(x['metabolite_key'],x['gene']) for x in old})==39
    interpretations={}
    for g in ['AQP9','SLC15A3','SLC29A3']:interpretations[g]='三研究宽髓系最高中位表达位置一致；不能推断同一髓系亚型、细胞特异性或代谢机制。'
    for g in ['BCAT2','UCKL1']:interpretations[g]='三研究上皮最高中位表达位置一致；上皮标签不等于确认恶性细胞。'
    for g in ['GSTA4','NNMT']:interpretations[g]='Lee/Pelka宽基质与Khaliq成纤维细胞线索相容；并非相同亚型复现，Khaliq仅7位成纤维细胞患者达标。'
    interpretations.update(SLC6A6='Lee髓系、Pelka基质、Khaliq成纤维细胞，保留多来源解释，不作多数票选择。',PRMT7='Lee与Khaliq上皮；Pelka v2为TNKILC、v3为上皮，保留技术/来源差异；基因级不能回答剪接变体。',KMT5A='Lee/Khaliq上皮与Pelka基质不一致，不能用多数票决定唯一来源。',HDC='Khaliq无肥大细胞标签，达标区室中位数均0；本批无法复核肥大细胞来源，不支持也不反证该来源。')
    for g in ['ASPA','ENTPD3','PGAM2','SLC38A3']:interpretations[g]='Khaliq达标区室中位CPM均0，不能据此定位；不等于所有细胞/患者均无表达，历史研究结果保留。'
    for g in ['GSTT2','SLC6A17','UPP2']:interpretations[g]='Khaliq作者区室导出缺独立条目，属当前矩阵不可评估；不同于低表达或零值，不扩大为全部原始资料不可测。'
    context={};genes=[]
    for gene in sorted({x['gene'] for x in old}):
        h=[x for x in historical if x['gene']==gene];k=[x for x in khaliq if x['gene']==gene];full=[x for x in primary if x['gene']==gene]
        assert len(h)==3 and len(k)==1 and len(full)==6
        assert {x['study'] for x in h}=={'Lee2020_GSE132465','Pelka2021_GSE178341'}
        views=[]
        for x in h:
            views.append({q:x[q] for q in ['study','technology','eligible_lineages','highest_median_CPM_lineage','n_qualified_patients','median_CPM','median_detection','status']})
        views.append(dict(study='Khaliq2022_GSE200997',technology='10x_5prime',**{q:k[0][q] for q in ['eligible_lineages','highest_median_CPM_lineage','n_qualified_patients','median_CPM','median_detection','status']}))
        detail=[{q:x[q] for q in ['cell_type','n','n_reference','n_cells','status','reason','pseudobulk_CPM_median','detection_fraction_median']} for x in full]
        item=dict(gene=gene,source3_version=RUN,source3_overview_json=json.dumps(views,ensure_ascii=False),source3_khaliq_coverage_expression_json=json.dumps(detail,ensure_ascii=False),source3_interpretation_cn=interpretations.get(gene,'保留四个技术展示层的原最高值和覆盖；不自动合并细胞标签或按多数票选择来源。'),source3_scope_cn='三项原始研究、四个技术展示层；Pelka v2/v3不是两项独立研究；仅已注释且覆盖达标类群的描述性最高值。',source3_metabolite_validation='NOT_RUN_BY_CELL_SOURCE_ANALYSIS',source3_khaliq_correction=CORRECTION)
        genes.append(item);context[gene]=item
    write(d/'gene35_source_context.tsv',list(genes[0]),genes)
    appended=[k for k in genes[0] if k!='gene'];assert not set(appended)&set(oldfields)
    new=[dict(x,**{k:context[x['gene']][k] for k in appended}) for x in old]
    write(d/'relations39_source_context.tsv',oldfields+appended,new)
    _,roundtrip=read(d/'relations39_source_context.tsv')
    assert len(roundtrip)==len(old)
    assert all(all(x[f]==y[f] for f in oldfields) for x,y in zip(old,roundtrip))
    admission=json.loads((r/IBASE/'admission_summary.json').read_text())
    assert admission['structural_pass']=={'RNA':30,'PROTEIN':11} and admission['correlations_computed']==0
    gates=[]
    for layer,n in [('RNA',30),('PROTEIN',11)]:
        gates.append(dict(layer=layer,structural_candidates=n,structural_status='DONE_REUSED',source_status='NEEDS_REVIEW',association_status='NOT_RUN',own_missing=('G_transcriptome_gene exact RSEM output;sample normalization;zero/imputation semantics' if layer=='RNA' else 'D_proteome_protein actual DIA-NN export field;sample normalization;NA/imputation semantics'),shared_missing='B_metabolome sample normalization and pre-export missing/imputation semantics',other_layer_blocks_this_layer=False,release_rule='Release this layer when its own and shared source questions are resolved;lock all eligible pairs before correlations',evidence=IBASE+'/README_CN.md'))
    write(d/'integral_layer_gates.tsv',list(gates[0]),gates)
    validation=dict(status='PASS',relations=39,genes=35,original_field_count=len(oldfields),all_original_fields_string_equal=True,original_row_order_preserved=True,original_relation_file_sha256=hashlib.sha256((r/OLD).read_bytes()).hexdigest(),studies=3,display_strata_per_gene=4,khaliq_compartments_per_gene=6,patient_statistics_recomputed=False,new_correlations=0,integral_structural_counts_reused=admission['structural_pass'])
    (d/'validation.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    spec=dict(analysis_version='COAD_cell_context3_v1',run_id=RUN,operation='append context columns only',statistical_unit='Not applicable: existing public aggregate integration',new_tests=False,family='NA',seed=None,source_comparison='Original labels retained;no subtype harmonization or majority vote',integral='RNA and protein independently gated;shared metabolome source gate retained',invariants=['original39 field strings and row order','974 candidates','9/21/5','all historical P/q'])
    (d/'analysis_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    manifest=[]
    for path in [OLD,LEE,KBASE+'/tumor_lineage_overview.tsv',KBASE+'/tumor_lineage_all35.tsv',IBASE+'/admission_summary.json',IBASE+'/README_CN.md']:
        manifest.append(dict(source_id=path,url_or_path=path,version='existing version preserved',sha256=hashlib.sha256((r/path).read_bytes()).hexdigest()))
    manifest += [dict(source_id='Khaliq_correction',url_or_path=CORRECTION,version='2022-07-13;Fig.S5 in Additional file1 only;accessed2026-09-21',sha256='NA:webpage not archived'),dict(source_id='Integral_methods_metadata',url_or_path='https://api.figshare.com/v2/articles/28156568',version='v1;accessed2026-09-21;no custom export fields',sha256='NA:web response not archived'),dict(source_id='Integral_S4_metadata',url_or_path='https://api.figshare.com/v2/articles/28156574',version='v1;accessed2026-09-21;no custom export fields',sha256='NA:web response not archived')]
    write(d/'source_manifest.tsv',list(manifest[0]),manifest)
    print(json.dumps(validation,ensure_ascii=False))

if __name__=='__main__':main()
