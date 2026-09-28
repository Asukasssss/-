"""Server-only donor-wise LYPLA1 coexpression and independent CRC context check.
Raw matrices and donor results never leave the server. No FDR selection.
"""
import argparse, gzip, io, json, tarfile, sys, hashlib, platform
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata, binomtest, spearmanr
ROOT=Path('/public3/xuzx/Cancer/pancancer_metabolomics_direct_matrix_20260716')
PRIOR=ROOT/'results/collaborative/COAD/B/20260925T140257Z_uckl1_cna_v1'
sys.path.insert(0,str(PRIOR))
from run_uckl1_cna_v1 import build_labels,sha
FOCUS='CDS1 ETNK2 GPAT4 GPD1L PISD PLA2G10 PLPP5 ASNS UCKL1 NNMT SLC6A6'.split()
PREFIX='cancer cohort stage_id run_id analysis_version analysis_type metabolite_key metabolite_name gene unit n n_reference effect_type effect ci_lower ci_upper p_value q_value test_family family_n_evaluable status reason source_id'.split()
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def write(p,x):pd.DataFrame(x).to_csv(p,sep='\t',index=False,na_rep='NA')
def ranknorm(a):
    z=rankdata(a,axis=0,method='average');return z-z.mean(axis=0)
def corr(x,y):
    num=x@y;den=np.sqrt((x*x).sum()*(y*y).sum(axis=0))
    return np.divide(num,den,out=np.full_like(num,np.nan,dtype=float),where=den>1e-10)
def residual(z,c):return z-c@np.linalg.lstsq(c,z,rcond=None)[0]
def signp(x):
    z=x[np.isfinite(x)];up=(z>1e-12).sum();down=(z< -1e-12).sum()
    return float(binomtest(int(up),int(up+down)).pvalue) if up+down else 1.
def normalized(c,t):return np.log1p(c/t[:,None]*10000)

def subtypes(meta,c,total,pub,out,study):
    x=np.log1p(c/total*10000);rows=[];private=[]
    for (donor,label),ids in meta.groupby(['donor','subtype']).groups.items():
        ids=np.array(list(ids));z=x[ids];det=c[ids]>0
        private.append(dict(donor=donor,subtype=label,cells=len(ids),mean=float(z.mean()),detection=float(det.mean()),positive_mean=float(z[det].mean()) if det.any() else np.nan))
    d=pd.DataFrame(private);write(out/('private_'+study+'_subtypes.tsv'),d)
    for label,ids in meta.groupby('subtype').groups.items():
        ids=np.array(list(ids));s=d[(d.subtype==label)&(d.cells>=20)];det=c[ids]>0
        rows.append(dict(cohort=study,subtype=label,cells=len(ids),donors_any=int(meta.loc[ids,'donor'].nunique()),donors_ge20=len(s),pooled_detection=float(det.mean()),pooled_mean=float(x[ids].mean()),pooled_positive_mean=float(x[ids][det].mean()) if det.any() else np.nan,donor_equal_mean=float(s['mean'].mean()) if len(s)>=3 else np.nan,donor_equal_detection=float(s.detection.mean()) if len(s)>=3 else np.nan,donor_equal_positive_mean=float(s.positive_mean.mean()) if len(s)>=3 else np.nan,status='DONE' if len(s)>=3 else 'NOT_EVALUABLE',reason='DESCRIPTIVE_AUTHOR_SUBTYPE' if len(s)>=3 else 'FEWER_THAN_3_DONORS_GE20'))
    write(pub/(study+'_subtype_summary.tsv'),rows)
    contrasts=[]
    for label in sorted(meta.subtype.unique()):
        dif=[]
        for donor,f in meta.groupby('donor'):
            a=f.index[f.subtype==label].to_numpy();b=f.index[f.subtype!=label].to_numpy()
            if min(len(a),len(b))>=20:dif.append(float(x[a].mean()-x[b].mean()))
        dx=np.array(dif);n=len(dx)
        contrasts.append(dict(subtype=label,contrast='subtype_vs_other_admitted_cells_within_donor',n_donors=n,median_difference=float(np.median(dx)) if n else np.nan,n_positive=int((dx>0).sum()),n_negative=int((dx<0).sum()),p_value=signp(dx) if n>=6 else np.nan,status='DONE' if n>=6 else 'NOT_EVALUABLE',reason='NOMINAL_EXPLORATORY' if n>=6 else 'FEWER_THAN_6_PAIRED_DONORS'))
    write(pub/(study+'_subtype_contrasts.tsv'),contrasts)
    return rows

def analyze(meta,C,genes,total,pub,out,study,spec):
    genes=np.asarray(genes);target=int(np.flatnonzero(genes=='LYPLA1')[0]);G=len(genes)
    donors=sorted(meta.donor.unique());modes=['raw','depth','depth_subtype','positive_only']
    matrices={mode:np.full((len(donors),G),np.nan) for mode in modes};coverage=[];checks=[]
    private=[]
    for di,donor in enumerate(donors):
        ix=meta.index[meta.donor==donor].to_numpy();counts=C[ix];t=total[ix];ly=counts[:,target]
        valid=len(ix)>=50 and (ly>0).sum()>=10 and np.unique(ly).size>=3
        coverage.append(dict(donor=donor,cells=len(ix),positive_cells=int((ly>0).sum()),included=bool(valid)))
        if not valid:continue
        x=ranknorm(np.log1p(ly/t*10000));depth=ranknorm(np.log1p(t-ly));controls=np.column_stack([np.ones(len(ix)),depth])
        # Author subtypes only, no expression-driven reclustering or label transfer.
        st=pd.get_dummies(meta.loc[ix,'subtype'],drop_first=True).to_numpy(dtype=float)
        cs=np.column_stack([controls,st]);pos=ly>0;positive_ok=pos.sum()>=30
        for start in range(0,G,256):
            stop=min(start+256,G);cc=counts[:,start:stop];v=normalized(cc,t);y=ranknorm(v)
            ok=((cc>0).sum(axis=0)>=10)&((y*y).sum(axis=0)>1e-10)
            for mode,xx,yy in [('raw',x,y),('depth',residual(x,controls),residual(y,controls)),('depth_subtype',residual(x,cs),residual(y,cs))]:
                r=corr(xx,yy);r[~ok]=np.nan;matrices[mode][di,start:stop]=r
            if positive_ok:
                xp=ranknorm(np.log1p(ly[pos]/t[pos]*10000));yp=ranknorm(v[pos]);rp=corr(xp,yp);rp[(cc[pos]>0).sum(axis=0)<10]=np.nan;matrices['positive_only'][di,start:stop]=rp
            if start==0:
                for k in np.flatnonzero(ok)[:5]:
                    expected=spearmanr(np.log1p(ly/t*10000),v[:,k]).statistic
                    checks.append(bool(np.isclose(matrices['raw'][di,k],expected,atol=1e-12)))
        print(study,'donor_done',di+1,'of',len(donors),flush=True)
    assert checks and all(checks)
    for mode,m in matrices.items():
        m[:,target]=np.nan
        for di,donor in enumerate(donors):
            private.extend(dict(donor=donor,gene=str(genes[g]),mode=mode,rho=float(m[di,g])) for g in np.flatnonzero(np.isfinite(m[di])))
    pd.DataFrame(private).to_csv(out/('private_'+study+'_correlations.tsv.gz'),sep='\t',index=False,compression='gzip')
    write(out/('private_'+study+'_coverage.tsv'),coverage)
    rows=[];rng=np.random.default_rng(spec['seed']);raw=matrices['raw']
    for g,gene in enumerate(genes):
        v=raw[:,g];v=v[np.isfinite(v)];n=len(v);row=dict.fromkeys(PREFIX,'NA')
        row.update(cancer='COAD',cohort=study,stage_id='06_EXTERNAL',run_id=out.name,analysis_version=spec['analysis_version'],analysis_type='donor_wise_LYPLA1_gene_rank_correlation',gene=gene,unit='donor',n=n,effect_type='median_within_donor_Spearman_rho',test_family=study+'_all_admitted_gene_correlations',source_id=study,status='DONE' if n>=6 else 'NOT_EVALUABLE',reason='NOMINAL_SIGN_TEST_ACROSS_DONORS' if n>=6 else 'FEWER_THAN_6_EVALUABLE_DONORS')
        if gene=='LYPLA1':row['reason']='SELF_CORRELATION_EXCLUDED'
        if n:
            row.update(effect=float(np.median(v)),n_positive=int((v>1e-12).sum()),n_negative=int((v< -1e-12).sum()),direction_fraction=float(max((v>1e-12).sum(),(v< -1e-12).sum())/n),rho_min=float(v.min()),rho_max=float(v.max()))
        if n>=6:
            boots=np.median(v[rng.integers(n,size=(2000,n))],axis=1);row.update(p_value=signp(v),ci_lower=float(np.quantile(boots,.025)),ci_upper=float(np.quantile(boots,.975)),leave_one_donor_out_min=float(min(np.median(np.delete(v,k)) for k in range(n))),leave_one_donor_out_max=float(max(np.median(np.delete(v,k)) for k in range(n))))
        for mode in modes[1:]:
            z=matrices[mode][:,g];z=z[np.isfinite(z)];row[mode+'_n']=len(z);row[mode+'_rho']=float(np.median(z)) if len(z) else np.nan;row[mode+'_p']=signp(z) if len(z)>=6 else np.nan
            row[mode+'_direction_fraction']=float(max((z>1e-12).sum(),(z< -1e-12).sum())/len(z)) if len(z) else np.nan
        rows.append(row)
    res=pd.DataFrame(rows);res['family_n_evaluable']=int((res.status=='DONE').sum())
    res['selected']=res.status.eq('DONE') & (pd.to_numeric(res.p_value,errors='coerce')<.05)&(pd.to_numeric(res.effect,errors='coerce').abs()>=.1)&(res.direction_fraction>=.75)
    res['depth_supported']=res.selected&(res.depth_n>=6)&(res.depth_p<.05)&(np.sign(pd.to_numeric(res.effect,errors='coerce'))==np.sign(res.depth_rho))
    write(pub/(study+'_results.tsv'),res);write(pub/(study+'_focus_genes.tsv'),res[res.gene.isin(FOCUS)])
    subtypes(meta,C[:,target],total,pub,out,study)
    dump(pub/(study+'_validation.json'),dict(status='PASS',cells=len(meta),donors=len(donors),donors_ge50_and_target_eligible=sum(x['included'] for x in coverage),genes=G,evaluable=int((res.status=='DONE').sum()),selected=int(res.selected.sum()),depth_supported=int(res.depth_supported.sum()),scipy_spearman_checks=len(checks),all_checks_passed=True,zero_values_retained=True,source_unit='cells nested within donor',no_fdr=True))
    return res

def extract_uhlitz(out,pub,spec):
    j,paths,_,audit=build_labels();assert {k:sha(p) for k,p in paths.items()}==json.loads((PRIOR/'admission_hashes.json').read_text())
    source=ROOT/'data/candidates/coad_uhlitz_20260921/counts.tar';paths['full_counts']=source
    expected=json.loads((ROOT/'results/collaborative/COAD/B/20260922T141755Z_source_contract_sc_v2/public/extraction_validation.json').read_text())['source_hashes']['counts.tar'];assert sha(source)==expected
    m=j[j.group=='tumor_CNA'].copy().reset_index(drop=True);m=m.rename(columns={'case_id':'donor','cell_type_epi_custom':'subtype'})
    assert len(m)==4477;m.to_csv(out/'private_Uhlitz_metadata.tsv',sep='\t',index=False)
    meta=pd.read_csv(paths['cell_metadata'],sep='\t',keep_default_na=False);ix=pd.Index(meta.cell_id).get_indexer(m.cell_id);assert (ix>=0).all()
    with np.load(paths['counts']) as z:total=z['total'][ix];ly=z['LYPLA1'][ix]
    genes=None;C=None;seen=set()
    with tarfile.open(source) as tar:
        for member in tar.getmembers():
            if not member.isfile():continue
            stream=io.TextIOWrapper(gzip.GzipFile(fileobj=tar.extractfile(member)));h=stream.readline().rstrip().split('\t');assert h[0]=='gene';order=pd.Index(m.cell_id).get_indexer(h[1:]);use=np.flatnonzero(order>=0)
            if not len(use):stream.close();continue
            cols=order[use];assert not seen.intersection(cols);seen.update(cols);names=[];values=[]
            for line in stream:
                name,raw=line.rstrip('\r\n').split('\t',1);v=np.fromstring(raw,sep='\t',dtype=np.int64)[use];assert (v>=0).all();names.append(name);values.append(v)
            stream.close();part=np.asarray(values,dtype=np.int32).T
            if genes is None:genes=names;C=np.zeros((len(m),len(genes)),dtype=np.int32)
            else:assert names==genes
            C[cols]=part;print('matrix_extracted',len(seen),'cells',flush=True)
    assert len(seen)==len(m) and len(set(genes))==len(genes);assert np.array_equal(C.sum(1),total);assert np.array_equal(C[:,genes.index('LYPLA1')],ly)
    np.savez_compressed(out/'private_Uhlitz_full_counts.npz',counts=C,genes=np.array(genes),total=total)
    dump(pub/'label_audit.json',audit);write(pub/'source_manifest.tsv',[dict(source_id=k,server_path=str(p),sha256=sha(p)) for k,p in paths.items()])
    return m,C,genes,total

def extract_lee(out,pub,targets):
    d=ROOT/'data/candidates/coad_cell_origin_20260920';annotation=d/'lee_annotation.txt.gz';source=d/'lee_raw_UMI.txt.gz'
    prev=ROOT/'results/collaborative/COAD/B/20260922T141755Z_source_contract_sc_v2/public/Lee_validation.json';hashes=json.loads(prev.read_text())['source_hashes']
    assert sha(annotation)==hashes[str(annotation)] and sha(source)==hashes[str(source)]
    a=pd.read_csv(annotation,sep='\t',dtype=str);mask=(a.Class=='Tumor')&(a.Cell_type=='Epithelial cells')&a.Cell_subtype.isin(['CMS1','CMS2','CMS3','CMS4'])
    m=a.loc[mask].rename(columns={'Index':'cell_id','Patient':'donor','Cell_subtype':'subtype'}).reset_index(drop=True);assert m.cell_id.is_unique
    targets=set(targets)|{'LYPLA1'};names=[];values=[];total=np.zeros(len(m),dtype=np.int64)
    with gzip.open(source,'rt') as f:
        h=[x.strip('"') for x in f.readline().rstrip().split('\t')];h=h[1:] if len(h)==len(a)+1 else h
        assert len(h)==len(a) and set(h)==set(a.Index);order=pd.Index(h).get_indexer(m.cell_id);assert (order>=0).all()
        for line in f:
            name,raw=line.rstrip('\r\n').split('\t',1);name=name.strip('"');v=np.fromstring(raw,sep='\t',dtype=np.int64)[order];assert (v>=0).all();total+=v
            if name in targets:names.append(name);values.append(v)
    assert len(names)==len(set(names)) and (total>0).all()
    old=ROOT/'results/collaborative/COAD/B/20260920T072000Z_cell_expression_v1/private_Lee_cell_UMI.npz'
    with np.load(old) as z:assert np.array_equal(total,z['total'][np.flatnonzero(mask)])
    C=np.array(values,dtype=np.int32).T;m.to_csv(out/'private_Lee_metadata.tsv',sep='\t',index=False)
    np.savez_compressed(out/'private_Lee_selected_counts.npz',counts=C,genes=np.array(names),total=total)
    write(pub/'Lee_feature_coverage.tsv',[dict(gene=g,status='DONE' if g in names else 'NOT_EVALUABLE',reason='EXACT_SYMBOL' if g in names else 'SOURCE_SYMBOL_ABSENT') for g in sorted(targets)])
    write(pub/'Lee_source_manifest.tsv',[dict(source_id=p.name,server_path=str(p),sha256=sha(p)) for p in [annotation,source,old]])
    return m,C,names,total

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();out=a.out;pub=out/'public';assert (out/'.running').is_dir();pub.mkdir(exist_ok=True)
    spec=json.loads((out/'analysis_spec.json').read_text());dump(pub/'analysis_spec.json',spec)
    m,C,genes,total=extract_uhlitz(out,pub,spec);res=analyze(m,C,genes,total,pub,out,'Uhlitz_CNA',spec)
    # Freeze replication targets using discovery only, before reading Lee expression.
    targets=sorted(set(res.loc[res.selected,'gene'])|set(FOCUS));dump(pub/'replication_targets_locked.json',dict(rule='all discovery selected + prespecified focus; no Lee expression used',genes=targets,n=len(targets)))
    del C
    m,C,genes,total=extract_lee(out,pub,targets);lee=analyze(m,C,genes,total,pub,out,'Lee_tumor_CMS_epithelium',spec)
    pairs=res[res.gene.isin(targets)].merge(lee,on='gene',how='outer',suffixes=('_Uhlitz','_Lee'),validate='one_to_one')
    pairs['direction_agrees']=np.sign(pd.to_numeric(pairs.effect_Uhlitz,errors='coerce'))==np.sign(pd.to_numeric(pairs.effect_Lee,errors='coerce'))
    pairs['nominal_replication']=pairs.selected_Uhlitz.fillna(False)&pairs.selected_Lee.fillna(False)&pairs.direction_agrees
    pairs['depth_supported_both']=pairs.nominal_replication & pairs.depth_supported_Uhlitz.fillna(False)&pairs.depth_supported_Lee.fillna(False)
    write(pub/'cross_CRC_comparison.tsv',pairs)
    dump(pub/'validation.json',dict(status='PASS',scope='Uhlitz donor-wise discovery and independent Lee tumor CMS epithelial context replication',replication_targets=len(targets),nominal_replication=int(pairs.nominal_replication.sum()),depth_supported_both=int(pairs.depth_supported_both.sum()),cross_cancer_status='NOT_RUN',Lee_does_not_replicate_CNA_call_or_TC_labels=True,software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__)))
    (out/'NUMERICAL_DONE').write_text('DONE\n');print('NUMERICAL_DONE',flush=True)
if __name__=='__main__':main()
