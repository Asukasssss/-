"""Small numeric fixture: normalization denominator, donor aggregation and missing genes."""
import argparse
import json
from pathlib import Path
import tempfile
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from brca_sc117_profile_v1 import run


def main(root):
    with tempfile.TemporaryDirectory(dir=root/'private',prefix='numeric_check_') as temp:
        p=Path(temp)
        for name in ['source','private','public']:(p/name).mkdir()
        genes=['G'+str(i) for i in range(117)]
        pd.DataFrame({'gene':genes}).to_csv(p/'source/candidate_comparison_v1.tsv',sep='\t',index=False)
        # Unequal cells/donor: equal-donor mean differs from a pooled cell mean.
        donors=[];labels=[];rows=[]
        for donor,m in [('D1',20),('D2',40),('D3',60)]:
            for celltype in ['A','B']:
                for j in range(m):
                    v=np.zeros(117,dtype=np.int64)
                    v[0]=int(donor[-1])*10 if celltype=='A' else 0
                    v[1]=20 if celltype=='A' else 50
                    v[-1]=100-v[0]-v[1]
                    rows.append(v);donors.append(donor);labels.append(celltype)
        x=sparse.csr_matrix(np.array(rows))
        with h5py.File(p/'source/test.h5ad','w') as h:
            dt=h5py.string_dtype('utf-8')
            o=h.create_group('obs');o.attrs['_index']='_index'
            for k,v in {'_index':['c'+str(i) for i in range(len(rows))], 'donor':donors,'label':labels}.items():
                o.create_dataset(k,data=np.array(v,dtype=object),dtype=dt)
            v=h.create_group('raw/var');v.attrs['_index']='_index'
            v.create_dataset('_index',data=np.array(genes[:-1]+['BG'],dtype=object),dtype=dt)
            v.create_dataset('feature_name',data=np.array(genes[:-1]+['BG'],dtype=object),dtype=dt)
            m=h.create_group('raw/X');m.attrs['shape']=x.shape;m.attrs['encoding-type']='csr_matrix'
            for k,a in [('data',x.data),('indices',x.indices),('indptr',x.indptr)]:m.create_dataset(k,data=a)
        config={'cohort':'numeric_fixture','file':'test.h5ad','raw_path':'raw/X','var_path':'raw/var',
                'symbol_column':'feature_name','donor_column':'donor','celltype_column':'label',
                'celltype_map':{'A':'A','B':'B'},'chunk_cells':17,'min_cells_per_donor_group':20,
                'min_donors_per_public_group':3,'unit':'synthetic_donor','source_id':'synthetic',
                'annotation_origin':'synthetic','independence_limit':'NA'}
        c=p/'config.json';c.write_text(json.dumps(config));run(p,c)
        d=pd.read_csv(p/'public/numeric_fixture_celltype_profiles.tsv',sep='\t')
        d=d[d.partition=='ALL']
        a=d[(d.gene=='G0')&(d.celltype=='A')].iloc[0]
        expected=float(np.mean(np.log1p(np.array([1000.,2000.,3000.]))))
        assert np.isclose(a.effect,expected,atol=1e-12)
        assert a['n']==3 and a.n_cells_total==120 and a.mean_detection_fraction==1
        b=d[(d.gene=='G0')&(d.celltype=='B')].iloc[0]
        assert b.effect==0 and b.mean_detection_fraction==0 and b.status=='DONE'
        absent=d[d.gene=='G116'];assert absent.effect.isna().all() and (absent.status=='NOT_EVALUABLE').all()
        assert d.gene.nunique()==117
        # Unambiguous identifiers are a gate, not silently de-duplicated.
        with h5py.File(p/'source/test.h5ad','r+') as h:h['obs/_index'][1]='c0'
        duplicate_rejected=False
        try:run(p,c)
        except AssertionError as e:duplicate_rejected='Nonunique cell' in str(e)
        assert duplicate_rejected
        report={'status':'PASS','unequal_donor_equal_weight_verified':True,'all_gene_library_denominator_verified':True,
                'zero_expression_distinguished_from_missing_gene':True,'all117_retained':True,
                'duplicate_barcodes_rejected':True,'expected_equal_donor_mean':expected,'observed':float(a.effect)}
        (root/'public/numeric_fixture_validation.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,required=True)
    main(parser.parse_args().root)
