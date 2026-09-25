"""Compact Chinese entry point for the complete 101-feature mapping ledger."""
import argparse,csv,json
from pathlib import Path


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);a=ap.parse_args();p=a.results
    names=json.loads(Path(__file__).with_name('paired_metabolite_names_cn.json').read_text(encoding='utf-8'))
    with (p/'feature_workpool101.tsv').open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    labels={'HAS_DIRECT_RELATIONS':'有直接注释关系','CONDITIONAL_IDENTITY':'身份/命名暂挂','CONDITIONAL_RELATIONS_ONLY':'仅条件关系','NO_DIRECT_HIT_IN_SCOPED_SEARCH':'本轮未找到直接支持'}
    md=['# COAD：101项配对代谢物的映射去向','',
        '78项形成891条主池关系、526基因；15项身份问题、2项仅条件关系、6项本轮无直接支持。主池是限定数据库注释支持的待检验集合，不是已验证基因名单。',
        '', '表中P/q为既有配对代谢物统计；新患者关联尚未计算。配对升/降/平的分母均33，作者填补限制保留。生成、消耗或转运角色见关系与证据表。','',
        '|代谢物|配对升/降/平|配对P|配对q|旧入口已有|去向|主池关系数|主池基因|',
        '|---|---|---:|---:|---|---|---:|---|']
    for r in rows:
        md.append('|{}|{}/{}/{}|{:.4g}|{:.4g}|{}|{}|{}|{}|'.format(names[r['metabolite_name']],r['up_pairs'],r['down_pairs'],r['equal_pairs'],
            float(r['paired_P']),float(r['paired_q']),'是' if r['legacy_feature']=='True' else '新增',labels[r['mapping_status']],r['direct_relations'],r['genes'] if r['genes'] not in ['', 'NA'] else '—'))
    (p/'FEATURE_OVERVIEW_CN.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
    print('Rendered101 feature rows')


if __name__=='__main__':main()
