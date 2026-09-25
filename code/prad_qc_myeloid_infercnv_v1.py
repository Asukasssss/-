"""Summarize selected-cell depth; individual measurements remain server-side."""
import argparse
from pathlib import Path
import pandas as pd,numpy as np
from scipy.io import mmread
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);r=p.parse_args().root
m=pd.read_csv(r/'private/selected_cell_metadata.tsv',sep='\t').set_index('cell');rows=[]
for donor in pd.read_csv(r/'private/donor_plan.tsv',sep='\t').donor_alias:
 d=r/'private'/donor;x=mmread(d/'counts.mtx').tocsc();cells=(d/'cells.txt').read_text().splitlines()
 counts=np.asarray(x.sum(axis=0)).ravel();genes=np.diff(x.indptr)
 rows.append(pd.DataFrame({'cell':cells,'mapped_autosomal_counts':counts,'detected_mapped_autosomal_genes':genes}).set_index('cell').join(m[['role']]))
z=pd.concat(rows);assert z.index.is_unique and set(z.index)==set(m.index)
z.to_csv(r/'private/input_cell_depth.tsv',sep='\t')
z['below100']=z.mapped_autosomal_counts.lt(100)
s=z.groupby('role').agg(n_cells=('mapped_autosomal_counts','size'),minimum_autosomal_counts=('mapped_autosomal_counts','min'),median_autosomal_counts=('mapped_autosomal_counts','median'),n_below100_counts=('below100','sum'),median_detected_autosomal_genes=('detected_mapped_autosomal_genes','median'))
s.to_csv(r/'public/06_EXTERNAL/input_depth_summary.tsv',sep='\t');print(s.to_string())
