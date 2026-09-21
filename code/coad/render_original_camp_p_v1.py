"""Chinese display of saved original COAD P values, not new tests."""
import argparse,csv,json,hashlib
from pathlib import Path


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);a=ap.parse_args();p=a.results
    names=json.loads(Path(__file__).with_name('paired_metabolite_names_cn.json').read_text(encoding='utf-8'))
    v=json.loads((p/'validation.json').read_text())
    for name,h in v['output_sha256'].items():assert hashlib.sha256((p/name).read_bytes()).hexdigest()==h
    with (p/'original_coad159.tsv').open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    chosen=sorted((r for r in rows if float(r['wilcoxon_p'])<.05),key=lambda r:float(r['wilcoxon_p']))
    assert len(chosen)==83
    md=['# COAD原CAMP代谢物：原P<0.05的全部83项','',
        '原159项中83项P<0.05，按原Hedges g方向分为40项肿瘤较高、43项肿瘤较低。73项同时原q<0.05，另10项仅名义P支持。',
        '', '这里直接读取原队列效应表wilcoxon_p，未重新计算P或q。原设计为37个肿瘤与39个正常标本的非配对比较，不是新做的33对患者检验，也不是代谢物—RNA关联。',
        '', 'g为标准化效应，不是倍数变化；正值表示肿瘤较高。＊标记原P<0.05但原q≥0.05的10项。中文名仅供阅读，英文原名及身份限制在TSV中保留。','']
    for positive,title in [(True,'肿瘤较高：40项'),(False,'肿瘤较低：43项')]:
        md+=['## '+title,'','|代谢物|原g|原P|原q|','|---|---:|---:|---:|']
        for r in chosen:
            if (float(r['hedges_g'])>0)==positive:
                md.append('|{}{}|{:+.3f}|{:.5g}|{:.5g}|'.format(names[r['feature_name']],'＊' if float(r['wilcoxon_fdr'])>=.05 else '',float(r['hedges_g']),float(r['wilcoxon_p']),float(r['wilcoxon_fdr'])))
        md.append('')
    (p/'ORIGINAL_P_LT_005_CN.md').write_text('\n'.join(md).rstrip()+'\n',encoding='utf-8')
    with (p/'original_P_lt_005_CN.tsv').open('w',encoding='utf-8',newline='') as f:
        fields=['代谢物','来源英文名','原g','原P','原q','原q是否小于0.05','肿瘤标本数','正常标本数','原分析设计']
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(fields)
        for r in chosen:w.writerow([names[r['feature_name']],r['feature_name'],r['hedges_g'],r['wilcoxon_p'],r['wilcoxon_fdr'],str(float(r['wilcoxon_fdr'])<.05),r['n_tumor'],r['n_normal'],r['analysis_design']])
    print('\n'.join(md))


if __name__=='__main__':main()
