"""Lossless cohort partitions for the repository's 5 MB file limit."""
import hashlib,json
import pandas as pd
from record_sop_v3_resume import SOURCE,write,write_json

def main():
    p=SOURCE/'sc_celltype_profiles.tsv';raw=p.read_bytes();lines=raw.splitlines(keepends=True);header=lines[0];cols=header.decode().strip().split('\t');idx=cols.index('cohort');groups={}
    for line in lines[1:]:groups.setdefault(line.decode().split('\t')[idx],[]).append(line)
    rows=[]
    for cohort,body in groups.items():
        f=SOURCE/f'sc_celltype_profiles_{cohort}.tsv';data=header+b''.join(body);assert len(data)<5_000_000;f.write_bytes(data)
        rows.append({'logical_table':'sc_celltype_profiles.tsv','cohort':cohort,'file':f.name,'rows':len(body),'sha256':hashlib.sha256(data).hexdigest(),'logical_table_sha256':hashlib.sha256(raw).hexdigest()})
    assert header+b''.join(line for body in groups.values() for line in body)==raw
    write(SOURCE/'sc_celltype_profile_parts.tsv',pd.DataFrame(rows));p.unlink()
    r=SOURCE/'README_CN.md';r.write_text(r.read_text()+'\n\n完整来源表为17,862行，按队列无损分为3个TSV以遵守单文件5 MB限制。字段、精度、顺序保持；分片清单及原逻辑表哈希见[sc_celltype_profile_parts.tsv](sc_celltype_profile_parts.tsv)。规范中的sc_celltype_profiles.tsv指这三个分片合并的逻辑表。\n',encoding='utf8',newline='\n')
    print(json.dumps({'lossless':True,'rows':sum(r['rows'] for r in rows),'parts':len(rows)}))
if __name__=='__main__':main()
