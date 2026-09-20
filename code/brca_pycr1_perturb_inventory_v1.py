"""Bounded source audit; ambiguous public metadata blocks formal inference."""
import argparse,hashlib,json,re
from pathlib import Path
import numpy as np
import pandas as pd

def main(root):
    src=root/'source';out=root/'public';out.mkdir(exist_ok=True)
    d=pd.read_csv(src/'GSE220931_processed_data.txt.gz',sep='\t')
    cols=list(d.columns[1:10]);a=d[cols].to_numpy(float)
    notes=[]
    for p in sorted(src.glob('GSM*.txt')):
        t=p.read_text();title=re.search(r'!Sample_title = (.+)',t).group(1).strip()
        genotype=re.search(r'genotype: (.+)',t).group(1).strip()
        notes.append(dict(accession=p.stem,title=title,genotype=genotype,in_matrix=title in cols))
    assert len(notes)==9 and all(x['in_matrix'] for x in notes)
    qc=dict(status='NEEDS_REVIEW',n_features=len(d),n_samples=len(cols),noninteger_values=int(np.count_nonzero(a!=np.floor(a))),missing_values=int(np.isnan(a).sum()),negative_values=int((a<0).sum()),duplicate_gene_ids=int(d.gene_id.duplicated().sum()),formal_tests=0,reason='Integer matrix without unambiguous count declaration; FPKM processing text; shP2_2/3 genotype labels conflict with titles; article hg19 vs GEO hg38',source_matrix_exported=False)
    (out/'validation.json').write_text(json.dumps(qc,indent=2)+'\n')
    pd.DataFrame(notes).to_csv(out/'sample_design_metadata.tsv',sep='\t',index=False)
    target=d.loc[d.gene_name.eq('PYCR1')];assert len(target)==1
    rows=[]
    for label in ['ctrl','shP1','shP2']:
        v=target[[c for c in cols if c.startswith(label+'_')]].to_numpy().ravel()
        rows.append(dict(group_as_column_title=label,n_columns=len(v),minimum_supplied_value=int(v.min()),median_supplied_value=float(np.median(v)),maximum_supplied_value=int(v.max()),unit='UNRESOLVED_AUTHOR_PROCESSED_UNITS',status='NEEDS_REVIEW',reason='Descriptive input audit only; not normalized abundance, knockdown efficiency or formal inference'))
    pd.DataFrame(rows).to_csv(out/'PYCR1_supplied_values_descriptive.tsv',sep='\t',index=False)
    pd.DataFrame([dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),server_path=str(p)) for p in sorted(src.iterdir()) if p.is_file()]).to_csv(out/'source_manifest.tsv',sep='\t',index=False)
    print(json.dumps(qc))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args().root)
