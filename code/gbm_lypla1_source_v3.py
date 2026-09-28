"""Server-only SNUH source QC and exact Cell Browser row extraction.

Run qc before extract; no target expression used for cohort selection.
"""
import argparse, hashlib, json, os, pathlib, struct, zlib
import numpy as np
import pandas as pd
import requests

BASE = 'https://cells.ucsc.edu/multiomic-gbm/scrna/'
MAP = {'tumor cell':'Malignant','microglia':'Myeloid','macrophage':'Myeloid',
       'monocyte':'Myeloid','DC':'Myeloid','granulocyte':'Myeloid',
       'lymphocyte':'Lymphoid','oligodendrocyte':'Oligodendrocyte',
       'endothelial':'Vascular','stromal cell':'Other_stromal'}
ALIASES = {'SNU21accutase':'SNU21','SNU33citeseq':'SNU33','SNU34citeseq':'SNU34'}

def save(d,p): d.to_csv(p,sep='\t',index=False,na_rep='NA')
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def row(root,gene):
    s=root/'source'; offsets=json.loads((s/'SNUH_exprMatrix.json').read_text())
    keys=[k for k in offsets if k==gene or k.split('|')[-1]==gene]
    assert len(keys)==1,(gene,keys)
    key=keys[0]; start,length=offsets[key][:2]; path=s/(gene+'.range.zlib')
    if not path.exists():
        proxy=os.environ.get('GBM_DOWNLOAD_PROXY')
        r=requests.get(BASE+'exprMatrix.bin?'+gene,headers={'Range':f'bytes={start}-{start+length-1}'},
                       proxies={'https':proxy} if proxy else None,timeout=(20,90))
        r.raise_for_status(); assert r.status_code==206
        assert r.headers['Content-Range'].startswith(f'bytes {start}-{start+length-1}/')
        assert len(r.content)==length;path.write_bytes(r.content)
    raw=zlib.decompress(path.read_bytes()); desc_len=struct.unpack('<H',raw[:2])[0]
    desc=raw[2:2+desc_len].decode(); values=np.frombuffer(raw[2+desc_len:],dtype='<f4').astype(float)
    assert len(values)==223113 and np.isfinite(values).all() and (values>=0).all()
    (s/(gene+'.range_provenance.json')).write_text(json.dumps(dict(url=BASE+'exprMatrix.bin',key=key,
        start=start,length=length,description=desc,n_cells=len(values),sha256=sha(path)),indent=2))
    return values

def metadata(root):
    d=pd.read_csv(root/'source/SNUH_meta_proxy.tsv',sep='\t',low_memory=False).rename(columns={'Unnamed: 0':'cell_id'})
    bar=pd.read_csv(root/'source/SNUH_barcodes.tsv.gz',header=None)[0]
    assert d.cell_id.is_unique and np.array_equal(d.cell_id.values,bar.values)
    d['patient']=d.ID1.replace(ALIASES)
    clinical=pd.read_excel(root/'source/SNUH_clinical_supplement.xlsx').set_index('patient_id')
    assert set(d.patient)<=set(clinical.index)
    expected=clinical[['scCore','scPeri','scNormal','csCore','csPeri','csNormal']].sum(axis=1)
    actual=d.groupby('patient').size()
    assert np.array_equal(actual.values,expected.loc[actual.index].values), 'clinical source counts disagree'
    d['category']=d['Major.celltype'].map(MAP); assert d.category.notna().all()
    keep=d.histology.eq('GBM') & d.Diagnosis.eq('Glioblastoma') & d.IDH1_mutation.eq('wildtype') & d.IDH2_mutation.eq('wildtype') & d.doblet.eq('Singlet') & d.tissue.isin(['core','peri','tumor'])
    assert (clinical.loc[d.loc[keep,'patient'],'Diagnosis']=='Glioblastoma').all()
    return d,keep

def main(root,mode):
    d,keep=metadata(root);s=root/'source';out=root/'public';out.mkdir(exist_ok=True)
    if mode=='qc':
        cohort=d[keep];summary=[]
        for cat,z in cohort.groupby('category'):
            counts=z.groupby('patient').size()
            summary.append(dict(category=cat,cells=len(z),patients=z.patient.nunique(),patients_ge20=int((counts>=20).sum()),
                patients_ge50=int((counts>=50).sum()),median_UMI=z.nCount_RNA.median(),median_genes=z.nFeature_RNA.median()))
        save(pd.DataFrame(summary),out/'SNUH_quality_before_expression.tsv')
        facts=dict(input_cells=len(d),gbm_before_singlet=int(d.histology.eq('GBM').sum()),selected_cells=int(keep.sum()),
                   patients=int(cohort.patient.nunique()),exact_barcode_order=True,clinical_patient_counts_exact=True,
                   author_aliases='three assay/preparation suffixes linked only after exact clinical count checks',
                   treatment_status='not supplied;do not label treatment-naive',
                   primary='all eligible whole-cell RNA,including CITE RNA;assay aliquots pooled by clinical patient',
                   sensitivities=['core_only','fresh_only','exclude_repeated_preparations','UMI1000_genes500'],
                   exclude_repeated_preparations='exclude SNU21accutase and SNU33citeseq;retain sole SNU34citeseq preparation',
                   normalization='author log(CPM+1);validate base with independent RPLP0 gene and metadata library totals')
        (s/'cohort_frozen_before_expression.json').write_text(json.dumps(facts,indent=2))
        print(json.dumps(facts));print(pd.DataFrame(summary).to_string(index=False))
        v=row(root,'RPLP0'); masks=(v>0)&(v<12);lib=d.nCount_RNA.values
        tests=[]
        for base in [np.e,2]:
            for factor in [1e4,1e6]:
                inv=(np.power(base,v[masks])-1)*lib[masks]/factor
                err=np.abs(inv-np.round(inv));tests.append(dict(base=float(base),factor=factor,median_integer_error=float(np.median(err)),fraction_within_001=float((err<.001).mean())))
        save(pd.DataFrame(tests),out/'normalization_probe.tsv');print(pd.DataFrame(tests).to_string(index=False))
        passed=[t for t in tests if t['base']==2 and t['factor']==1e6][0]
        assert passed['fraction_within_001']>.99
        (s/'normalization_verified.json').write_text(json.dumps(dict(scale='author log2(CPM+1)',
            evidence='paper CPM+log; independent RPLP0 inverse transform matches integer counts using nCount_RNA in >99% positive cells',
            no_additional_transform=True,**passed),indent=2))
    else:
        assert (s/'cohort_frozen_before_expression.json').exists() and (s/'normalization_verified.json').exists()
        v=row(root,'LYPLA1'); d['expression']=v;d['detected']=(v>0).astype(int)
        save(d.loc[keep,['cell_id','patient','category','Major.celltype','celltype','tissue','source','ID1','nCount_RNA','nFeature_RNA','expression','detected']],root/'private/SNUH_LYPLA1_cells.tsv')
        profiles=[]
        for partition,m in [('all',keep),('core_only',keep&d.tissue.eq('core')),('fresh_only',keep&d.source.eq('Fresh')),
            ('exclude_repeated_preparations',keep&~d.ID1.isin(['SNU21accutase','SNU33citeseq'])),
            ('UMI1000_genes500',keep&d.nCount_RNA.ge(1000)&d.nFeature_RNA.ge(500))]:
            z=d[m];a=z.groupby(['patient','category']).agg(n_cells=('expression','size'),sum_expression=('expression','sum'),sum_detected=('detected','sum')).reset_index()
            a['cohort']='SNUH';a['partition']=partition;profiles.append(a)
        save(pd.concat(profiles,ignore_index=True),root/'private/SNUH_patient_profiles.tsv')
        print('Target extraction finished; patient/cell values retained server-side.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=pathlib.Path);p.add_argument('mode',choices=['qc','extract']);a=p.parse_args();main(a.root,a.mode)
