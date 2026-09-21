"""Display existing paired P<0.05 by direction counts; no new statistical tests."""
import argparse,csv,hashlib,json
from collections import defaultdict
from pathlib import Path


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=False)
    namespath=Path(__file__).with_name('paired_metabolite_names_cn.json');names=json.loads(namespath.read_text(encoding='utf-8'))
    v=json.loads((a.source/'validation.json').read_text())
    source=a.source/'primary33.tsv';assert sha(source)==v['output_sha256']['primary33.tsv']
    with source.open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==159 and all(int(r['n'])==33 and sum(int(r[k]) for k in ['up_pairs','down_pairs','equal_pairs'])==33 for r in rows)
    chosen=[r for r in rows if float(r['p_value'])<.05]
    summary=dict(all_features=159,paired_P_lt_005=len(chosen),majority_up=sum(int(r['up_pairs'])>16 for r in chosen),
        majority_down=sum(int(r['down_pairs'])>16 for r in chosen),no_majority=sum(max(int(r['up_pairs']),int(r['down_pairs']))<=16 for r in chosen),
        q_lt_005=sum(float(r['q_value'])<.05 for r in chosen),nominal_only=sum(float(r['q_value'])>=.05 for r in chosen),new_tests=False)
    assert summary['paired_P_lt_005']==101 and summary['nominal_only']==7
    md=['# COAD：33对患者中配对P<0.05的结果，按同向配对率排列','',
        '共101项：46项多数升高、54项多数降低、1项没有过半患者同向。全部沿用既有33对Wilcoxon配对P，不使用原37/39非配对P，未重算P/q。',
        '', '比例=升高或降低对数÷33，相等对也保留在分母内。＊表示配对P<0.05但配对q≥0.05；†表示至少一对涉及作者填补值。中文名仅作显示，来源身份限制保留。','']
    for field,title in [('up_pairs','多数升高：46项'),('down_pairs','多数降低：54项')]:
        grouped=defaultdict(list)
        for r in chosen:
            if int(r[field])>16:
                grouped[int(r[field])].append(names[r['metabolite_name']]+('＊' if float(r['q_value'])>=.05 else '')+('†' if int(r['both_observed_pairs'])<33 else ''))
        md+=['## '+title,'','|同向对数|配对率|代谢物——每项均达到该行对数|','|---:|---:|---|']
        for n,labels in sorted(grouped.items(),reverse=True):md.append('|{}/33|{:.1f}%|{}|'.format(n,n/33*100,'；'.join(labels)))
        md.append('')
    md+=['## 单列：葡萄糖醛酸','',
        '7/33升高（21.2%）、16/33降低（48.5%）、10/33相等（30.3%）；配对P=0.026399，q=0.044654。不能称多数降低。','',
        '## 7项仅配对P显著、q未显著','', '|代谢物|升高/降低/相等|配对P|配对q|','|---|---|---:|---:|']
    for r in chosen:
        if float(r['q_value'])>=.05:md.append('|{}|{}/{}/{}|{:.6g}|{:.6g}|'.format(names[r['metabolite_name']],r['up_pairs'],r['down_pairs'],r['equal_pairs'],float(r['p_value']),float(r['q_value'])))
    (a.out/'PAIRED_P_COUNTS_CN.md').write_text('\n'.join(md).rstrip()+'\n',encoding='utf-8')
    # Preserve all numerical strings; change only the delivery identifiers and add display fields.
    fields=list(rows[0])+['source_run_id','source_analysis_version','metabolite_name_cn','paired_P_lt_005','paired_q_lt_005']
    with (a.out/'results.tsv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for r in rows:
            x=dict(r,source_run_id=r['run_id'],source_analysis_version=r['analysis_version'],run_id=a.out.name,
                analysis_version='paired_p_counts_view_v1',metabolite_name_cn=names[r['metabolite_name']],
                paired_P_lt_005=str(float(r['p_value'])<.05),paired_q_lt_005=str(float(r['q_value'])<.05))
            assert all(x[k]==r[k] for k in r if k not in ['run_id','analysis_version'])
            w.writerow(x)
    for file,obj in [('summary.json',summary),('analysis_spec.json',dict(version='paired_p_counts_view_v1',source_run=a.source.name,selection='existing primary33 p_value < 0.05',ordering='up/down count descending, denominator 33 incl equal',new_tests=False)),
        ('validation.json',dict(status='PASS',source_hash_matches=True,all159_count_sums33=True,all_source_fields_preserved_except_view_identifiers=True,selected101=True,all101_accounted_for=True))]:
        (a.out/file).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    with (a.out/'source_manifest.tsv').open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['source_path','version','sha256'])
        for p in [source,a.source/'validation.json',namespath,Path(__file__)]:w.writerow([p.as_posix(),'existing_aggregate_or_display_code',sha(p)])
    print('\n'.join(md))


if __name__=='__main__':main()
