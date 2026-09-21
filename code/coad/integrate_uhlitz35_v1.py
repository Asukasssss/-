"""Append fourth-study aggregates; preserve every historical relation field."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

OLD='results/COAD/07_INTEGRATION/20260921T095621Z_cell_context3_v1/relations39_source_context.tsv'

def read(p):
    with p.open(encoding='utf-8',newline='') as f:
        r=csv.DictReader(f,delimiter='\t');return r.fieldnames,list(r)

def write(p,fields,rows):
    with p.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--batch',type=Path,required=True);a=p.parse_args()
    d=a.batch;fields,old=read(a.repo/OLD);_,summary=read(d/'tumor_overview.tsv');_,detail=read(d/'tumor_all35.tsv')
    assert len(old)==39 and len(set(x['gene'] for x in old))==35
    context=[]
    notes={
        'HDC':'本批肥大细胞8位患者达标，新增肥大细胞表达支持；Khaliq缺少该注释仍不构成反证，不验证组胺关联。',
        'GSTA4':'本批周细胞最高、内皮次之，宽基质背景相容；不能概括为四研究成纤维细胞特异来源。',
        'BCAT2':'本批周细胞中位CPM略高于上皮；周细胞5位、上皮12位患者达标，非相同患者集合的配对比较。不能写四研究上皮最高一致。',
        'NNMT':'本批成纤维细胞10位患者达标并为最高中位CPM类群，增加基质/成纤维背景支持；非唯一来源。',
        'AQP9':'本批作者Monocyte类群12位患者达标且中位表达最高，与历史宽髓系线索相容；Monocyte不等于全髓系。',
        'UCKL1':'本批上皮12位患者达标且中位表达最高，与此前上皮线索相容；不自动等于恶性细胞特异性。',
        'SLC6A6':'本批内皮中位表达最高，保留免疫与基质多来源解释；不强行归一来源。',
        'PRMT7':'本批上皮最高，但保留历史研究技术/来源差异；基因表达仍不能定量PRMT7变体。',
        'KMT5A':'本批上皮与内皮中位CPM接近，不能仅凭最高标签排除内皮或历史基质背景。'}
    for gene in ['GSTT2','GSTT2B','SLC6A17','UPP2']:
        notes[gene]='本批作者GEO导出矩阵缺少可用独立条目，不是零表达；不扩大解释为原始测序或该队列无法测量。'
    notes['SLC38A3']='本批上皮最高，但患者检出率中位数仅约0.0939%，仍属低检出线索，不能仅靠最高标签宣布明确来源。'
    notes['SLC29A3']='本批Monocyte最高，与历史宽髓系线索相容，但检出率中位数仅约0.7083%，须保留低检出限制。'
    for gene in sorted(set(x['gene'] for x in old)):
        overview=[x for x in summary if x['gene']==gene];rows=[x for x in detail if x['gene']==gene]
        assert len(overview)==3
        keep=['annotation_level','cell_type','n','n_reference','n_cells','status','pseudobulk_CPM_median',
              'pseudobulk_CPM_q25','pseudobulk_CPM_q75','detection_fraction_median']
        context.append(dict(gene=gene,source4_run=d.name,source4_uhlitz_overview_json=json.dumps(overview,ensure_ascii=False),
            source4_uhlitz_expression_json=json.dumps([{k:x[k] for k in keep} for x in rows],ensure_ascii=False),
            source4_interpretation_cn=notes.get(gene,'保留各原注释和预定合并类群的实际分布与覆盖；不以最高值自动宣布来源一致、特异或改变功能分层。'),
            source4_scope='Fourth original recruitment study;derived compartments plus original main/minor labels;no majority-vote source assignment',
            source4_boundary='Cell expression support only;no new metabolite association,paired differential test,P/q or change to9/21/5'))
    write(d/'gene35_source4_context.tsv',list(context[0]),context)
    by={x['gene']:x for x in context};extra=[x for x in context[0] if x!='gene'];assert not(set(extra)&set(fields))
    new=[dict(x,**{k:by[x['gene']][k] for k in extra}) for x in old]
    write(d/'relations39_source4_context.tsv',fields+extra,new)
    _,reread=read(d/'relations39_source4_context.tsv')
    assert all(all(x[k]==y[k] for k in fields) for x,y in zip(old,reread)) and len(reread)==39
    result=dict(status='PASS',relations=39,genes=35,original_fields=len(fields),original_order_and_strings_unchanged=True,
        source_path=OLD,source_sha256=hashlib.sha256((a.repo/OLD).read_bytes()).hexdigest(),new_p_q=False)
    (d/'integration_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))

if __name__=='__main__':main()
