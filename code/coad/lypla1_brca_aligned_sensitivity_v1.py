"""COAD-only numerical sensitivity aligned to existing BRCA definitions."""
import argparse,gzip,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import binomtest
from lypla1_state_coexpression_v1 import ROOT,sha
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();o=a.out;p=o/'public';assert (o/'NUMERICAL_DONE').exists()
    targets=json.loads((p/'replication_targets_locked.json').read_text())['genes'];coverage=[];allrows=[]
    for co,filename in [('Uhlitz','private_Uhlitz_full_counts.npz'),('Lee','private_Lee_selected_counts.npz')]:
        m=pd.read_csv(o/('private_'+co+'_metadata.tsv'),sep='\t',dtype={'donor':str},keep_default_na=False)
        with np.load(o/filename) as z:C=z['counts'];genes=z['genes'];total=z['total']
        target=int(np.flatnonzero(genes=='LYPLA1')[0]);ly=C[:,target].copy();keep=np.flatnonzero(np.isin(genes,targets));gene_names=genes[keep]
        if co=='Uhlitz':nfeat=(C>0).sum(1)
        else:
            source=ROOT/'data/candidates/coad_cell_origin_20260920/lee_raw_UMI.txt.gz';nfeat=np.zeros(len(m),int);den=np.zeros(len(m),np.int64)
            with gzip.open(source,'rt') as f:
                h=[x.strip('"') for x in f.readline().rstrip().split('\t')];h=h[1:] if h[0] not in set(m.cell_id) and len(h)==63690 else h
                ix=pd.Index(h).get_indexer(m.cell_id);assert (ix>=0).all()
                for line in f:
                    _,v=line.rstrip('\r\n').split('\t',1);v=np.fromstring(v,sep='\t',dtype=np.int64)[ix];nfeat+=v>0;den+=v
            assert np.array_equal(den,total)
        X=C[:,keep];del C
        t=np.log1p(ly/total*10000);depth=np.log1p(total-ly);det=np.log1p(nfeat-(ly>0));donors=[]
        for donor,ix in m.groupby('donor').groups.items():
            ix=np.array(list(ix));valid=len(ix)>=100 and (ly[ix]>0).sum()>=10 and np.std(t[ix])>0
            if valid:donors.append((donor,ix))
        admitted=np.concatenate([x[1] for x in donors]);globalok=(X[admitted]>0).mean(0)>=.01
        A={mode:np.full((len(donors),len(keep)),np.nan) for mode in ['raw','depth','depth_state']}
        for di,(donor,ix) in enumerate(donors):
            y=np.log1p(X[ix]/total[ix,None]*10000);x=t[ix];q=np.column_stack([np.ones(len(ix)),depth[ix],det[ix]])
            states=pd.get_dummies(m.loc[ix,'subtype'],drop_first=True).to_numpy(float);qs=np.column_stack([q,states])
            ok=globalok & ((X[ix]>0).sum(0)>=max(10,.01*len(ix)))
            for mode,c in [('raw',q[:,:1]),('depth',q),('depth_state',qs)]:
                rx=x-c@np.linalg.lstsq(c,x,rcond=None)[0];ry=y-c@np.linalg.lstsq(c,y,rcond=None)[0]
                vx=(rx*rx).sum();vy=(ry*ry).sum(0);good=ok&(vy>1e-10)&(vx>1e-10)
                A[mode][di,good]=np.clip((rx@ry[:,good])/np.sqrt(vx*vy[good]),-1,1)
            print('aligned',co,'donor',di+1,'/',len(donors),flush=True)
        for mode,z in A.items():
            pd.DataFrame(z,index=[d for d,_ in donors],columns=gene_names).to_csv(o/('private_'+co+'_aligned_'+mode+'.tsv.gz'),sep='\t',compression='gzip')
            for gi,g in enumerate(gene_names):
                v=z[:,gi];v=v[np.isfinite(v)];n=len(v);up=int((v>0).sum());down=int((v<0).sum())
                allrows.append(dict(cancer='COAD',cohort=co,mode=mode,gene=g,n_patients=n,median_r=float(np.median(v)) if n else np.nan,positive_fraction=up/n if n else np.nan,negative_fraction=down/n if n else np.nan,p_value=float(binomtest(up,up+down).pvalue) if n>=5 and up+down else np.nan,status='DONE' if n>=5 else 'NOT_EVALUABLE'))
        coverage.append(dict(cohort=co,cells_total=len(m),cells_eligible=len(admitted),donors_eligible=len(donors),genes_requested=len(targets),genes_present=len(keep),gene_detection_filter=.01))
    pd.DataFrame(allrows).to_csv(p/'COAD_BRCA_aligned_correlations.tsv',sep='\t',index=False,na_rep='NA')
    (p/'aligned_sensitivity_spec.json').write_text(json.dumps(dict(method='Pearson of log1p(CP10K); intercept + log1p(total UMI minus LYPLA1) + log1p(nDetected minus LYPLA1 detection); state dummies sensitivity',admission='>=100 cells,>=10 LYPLA1 positive,nonconstant target; gene>=1% eligible cells in cohort and >=max(10,1% cells) per donor;>=5 donors',scope='Frozen Uhlitz discovery targets plus focus, not all transcriptome; no BRCA matrices rerun',selection='No P or FDR gate for cross-cancer descriptive comparison; report >=0.1 effect and >=70% direction separately',code_sha256=sha(Path(__file__)),coverage=coverage),indent=2)+'\n')
    (o/'ALIGNED_DONE').write_text('DONE\n')
if __name__=='__main__':main()
