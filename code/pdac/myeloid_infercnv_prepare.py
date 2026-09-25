"""Extract all-gene inputs on server; never exports cells or donor profiles."""
from pathlib import Path
import gzip,hashlib,json
import numpy as np,pandas as pd
from scipy import io,sparse
RUN=Path(__file__).resolve().parent
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
DATA=ROOT/'data/candidates'
ORDER=pd.read_csv(RUN/'private/gene_order_no_chr3.tsv',sep='\t',header=None,names=['gene','chr','start','end'])
SEED=20260925
def seed(x):return int(hashlib.sha256(x.encode()).hexdigest()[:8],16)
def lines(p):
    with gzip.open(p,'rt') as f:return [x.rstrip('\n\r') for x in f]
def selection(c):
    d=pd.read_csv(RUN/'private'/f'{c}_metadata.tsv',sep='\t',dtype={'unit':str})
    counts=pd.read_csv(RUN/'private'/f'{c}_unit_counts.tsv',sep='\t',dtype={'unit':str})
    units=counts.loc[counts.eligible,'unit'];out=[]
    for unit in units:
        for role,g in d[d.unit==unit].groupby('role'):
            g=g.sort_values('cell').copy();rng=np.random.default_rng(seed(c+unit+role))
            if role in ['T_ref','B_ref']:
                idx=rng.permutation(len(g))[:200];g=g.iloc[idx].copy();n=len(g)//2
                g['cnv_label']=[role]*n+[role+'_heldout']*(len(g)-n)
            elif role in ['Author_malignant_control','Myeloid_normal_control']:
                g=g.iloc[rng.permutation(len(g))[:100]].copy();g['cnv_label']=role
            else:g['cnv_label']='Myeloid'
            out.append(g)
    return pd.concat(out).set_index('cell') if out else d.iloc[:0].set_index('cell')
def emit(c,meta,x,features):
    assert x.shape==(len(features),len(meta)) and np.all(x.data>=0)
    features=pd.Index(features);unique=~features.duplicated(keep=False)
    si=np.flatnonzero((features=='SLC6A6')&unique);assert len(si)==1
    lib=np.asarray(x.sum(axis=0)).ravel();assert (lib>0).all()
    slc=x[si[0]].toarray().ravel();meta=meta.copy()
    meta['slc6a6_log1p10k']=np.log1p(10000*slc/lib);meta['slc6a6_detected']=(slc>0).astype(int)
    meta['total_umi']=lib;meta['n_genes']=np.asarray((x>0).sum(axis=0)).ravel()
    # Author-QC retained; additional low-complexity filter fixed before SLC6A6 comparisons.
    meta['qc_pass']=(meta.total_umi>=500)&(meta.n_genes>=200)
    ix={g:i for i,g in enumerate(features) if unique[i]};genes=[g for g in ORDER.gene if g in ix]
    geneidx=np.array([ix[g] for g in genes]);assert 'SLC6A6' not in genes
    records=[]
    for unit,g in meta.groupby('unit',sort=True):
        tag=hashlib.sha256((c+unit).encode()).hexdigest()[:12]
        dest=RUN/'private/units'/c/tag;dest.mkdir(parents=True,exist_ok=True)
        g.to_csv(dest/'all_selected_metadata.tsv',sep='\t',index_label='cell')
        g=g[g.qc_pass];counts=g.cnv_label.value_counts()
        ok=counts.get('Myeloid',0)>=50 and all(counts.get(k,0)>=20 for k in ['T_ref','B_ref','T_ref_heldout','B_ref_heldout'])
        (dest/'unit.json').write_text(json.dumps({'cohort':c,'unit':str(unit),'tag':tag,'status':'READY' if ok else 'NOT_EVALUABLE','reason':'PostQC at least50 myeloid and20 each train/heldout T/B'},indent=2))
        if ok:
            pos=meta.index.get_indexer(g.index);assert (pos>=0).all()
            io.mmwrite(str(dest/'counts.mtx'),x[geneidx][:,pos].tocoo(),field='integer')
            (dest/'genes.txt').write_text('\n'.join(genes)+'\n')
            (dest/'cells.txt').write_text('\n'.join(g.index)+'\n')
            g.to_csv(dest/'metadata.tsv',sep='\t',index_label='cell')
            g[['cnv_label']].to_csv(dest/'annotations.tsv',sep='\t',header=False)
        records.append({'cohort':c,'unit':unit,'path':str(dest),'ready':ok,'myeloid':int(counts.get('Myeloid',0))})
    return records
def run242(meta):
    folder=DATA/'GSE242230/raw_matrices';records=[]
    for unit,g in meta.groupby('unit',sort=True):
        pp=list(folder.glob('*_'+unit+'_matrix.mtx.gz'));assert len(pp)==1
        p=pp[0];pre=str(p)[:-len('matrix.mtx.gz')]
        bars=[unit+'_'+v for v in lines(Path(pre+'barcodes.tsv.gz'))]
        features=[v.split('\t')[1] for v in lines(Path(pre+'features.tsv.gz'))]
        pos=pd.Index(bars).get_indexer(g.index);assert (pos>=0).all()
        with gzip.open(p,'rb') as f:x=io.mmread(f).tocsr()[:,pos]
        records+=emit('GSE242230',g,x,features);print('prepared GSE242230 unit',len(records),flush=True)
    return records
def run278(meta):
    folder=DATA/'gse278688_pdac_scRNA_v0.1'
    bars=lines(folder/'GSE278688_sc_barcodes.tsv.gz');features=[v.split('\t')[0] for v in lines(folder/'GSE278688_sc_features.tsv.gz')]
    pos=pd.Index(bars).get_indexer(meta.index);assert (pos>=0).all()
    cmap=np.full(len(bars),-1,dtype=np.int32);cmap[pos]=np.arange(len(pos));rr=[];cc=[];vv=[]
    with gzip.open(folder/'GSE278688_sc_matrix.mtx.gz','rt') as f:
        assert f.readline().startswith('%%MatrixMarket');line=f.readline()
        while line.startswith('%'):line=f.readline()
        nr,nc,nnz=map(int,line.split());assert(nr,nc)==(len(features),len(bars));seen=0
        for b in pd.read_csv(f,sep=r'\s+',header=None,names=['r','c','v'],dtype=np.int64,chunksize=2000000):
            a=b.to_numpy();seen+=len(a);ci=cmap[a[:,1]-1];keep=ci>=0
            rr.append((a[keep,0]-1).astype(np.int32));cc.append(ci[keep]);vv.append(a[keep,2].astype(np.int32))
        assert seen==nnz
    x=sparse.coo_matrix((np.concatenate(vv),(np.concatenate(rr),np.concatenate(cc))),shape=(len(features),len(meta))).tocsr()
    return emit('GSE278688',meta,x,features)
def run263(meta):
    folder=DATA/'gse263733_pdac_scRNA_v0.1';rr=[];cc=[];vv=[];genes=[]
    with gzip.open(folder/'GSE263733_Raw_counts.txt.gz','rt') as f:
        bars=f.readline().rstrip().split('\t');pos=pd.Index(bars).get_indexer(meta.index);assert (pos>=0).all()
        for line in f:
            name,values=line.split('\t',1);a=np.fromstring(values,sep='\t',dtype=np.int64);assert len(a)==len(bars)
            a=a[pos];nz=np.flatnonzero(a);rr.append(np.full(len(nz),len(genes),dtype=np.int32));cc.append(nz.astype(np.int32));vv.append(a[nz].astype(np.int32));genes.append(name)
    x=sparse.coo_matrix((np.concatenate(vv),(np.concatenate(rr),np.concatenate(cc))),shape=(len(genes),len(meta))).tocsr()
    return emit('GSE263733',meta,x,genes)
if __name__=='__main__':
    allrecords=[]
    for c,fn in [('GSE242230',run242),('GSE278688',run278),('GSE263733',run263)]:
        prior=RUN/'private/units_index.tsv'
        cached=pd.read_csv(prior,sep='\t') if prior.exists() else pd.DataFrame()
        if c=='GSE242230' and len(cached) and (cached.cohort==c).sum()==9:
            allrecords+=cached[cached.cohort==c].to_dict('records')
        else:allrecords+=fn(selection(c))
        pd.DataFrame(allrecords).to_csv(RUN/'private/units_index.tsv',sep='\t',index=False)
        print('COHORT_PREPARED',c,flush=True)
    d=pd.DataFrame(allrecords);summary=d.groupby('cohort').agg(n_units=('ready','size'),n_ready=('ready','sum'),n_myeloid_postQC=('myeloid','sum')).reset_index()
    summary.to_csv(RUN/'public/prepared_summary.tsv',sep='\t',index=False)
    (RUN/'PREPARED.json').write_text(summary.to_json(orient='records'))
