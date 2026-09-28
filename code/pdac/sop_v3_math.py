"""Small, explicitly defined SOP statistics; no I/O or gene selection."""
import hashlib
import numpy as np

def seed_for(cohort,family,key,master=20260922):
    return int(hashlib.sha256(f'{master}|PDAC_CAMP_source_SOP_v3|{cohort}|{family}|{key}'.encode()).hexdigest()[:8],16)

def bh_evaluable(values):
    p=np.asarray(values,float);valid=np.isfinite(p);q=np.full(len(p),np.nan)
    if np.any((p[valid]<0)|(p[valid]>1)):raise ValueError('P outside[0,1]')
    if valid.any():
        pv=p[valid];order=np.argsort(pv);ranked=pv[order]*len(pv)/np.arange(1,len(pv)+1)
        z=np.empty(len(pv));z[order]=np.minimum(1,np.minimum.accumulate(ranked[::-1])[::-1]);q[valid]=z
    return q

def paired_t(delta,seed,min_pairs=8,bootstrap_B=4000):
    from scipy.stats import t as student_t
    d=np.asarray(delta,float);d=d[np.isfinite(d)];n=len(d)
    z={'n':n,'n_up':int((d>0).sum()),'n_down':int((d<0).sum()),'n_equal':int((d==0).sum()),'mean_delta':float(d.mean()) if n else np.nan,'median_delta':float(np.median(d)) if n else np.nan,'p_value':np.nan,'ci_lower':np.nan,'ci_upper':np.nan,'bootstrap_mean_lower':np.nan,'bootstrap_mean_upper':np.nan,'bootstrap_valid':0,'status':'NOT_EVALUABLE','reason':'N_PAIRS_LT_8'}
    if n<min_pairs:return z
    if np.var(d,ddof=1)==0:z['reason']='ZERO_VARIANCE_PAIRED_DIFFERENCES';return z
    se=d.std(ddof=1)/np.sqrt(n);stat=d.mean()/se;crit=student_t.ppf(.975,n-1)
    rng=np.random.default_rng(seed);boot=d[rng.integers(n,size=(bootstrap_B,n))].mean(axis=1)
    z.update(status='DONE',reason='Paired t on author continuous scale',p_value=float(2*student_t.sf(abs(stat),n-1)),ci_lower=float(d.mean()-crit*se),ci_upper=float(d.mean()+crit*se),bootstrap_mean_lower=float(np.quantile(boot,.025)),bootstrap_mean_upper=float(np.quantile(boot,.975)),bootstrap_valid=bootstrap_B,t_statistic=float(stat),df=n-1)
    return z

def cellwise_profiles(counts,library,groups,present):
    """Rows genes, columns cells. groups contains cell index arrays per donor/type."""
    x=np.asarray(counts);lib=np.asarray(library,float)
    if x.ndim!=2 or x.shape[1]!=len(lib) or (x<0).any() or (lib<=0).any() or not np.isfinite(lib).all():raise ValueError('Invalid counts/library')
    if np.any(x.max(axis=0)>lib):raise ValueError('Gene count exceeds whole-library count')
    z=np.log1p(x*10000/lib[None,:]);out=[]
    for idx in groups:
        idx=np.asarray(idx,int);mean=z[:,idx].mean(axis=1);det=(x[:,idx]>0).mean(axis=1)
        mean[~present]=np.nan;det[~present]=np.nan
        out.append((len(idx),mean,det))
    return out

def source_rank_bootstrap(expression,detection,seed,B=1000):
    """Donor rows, type columns; NaN excluded donor/type, joint donor bootstrap."""
    x=np.asarray(expression,float);det=np.asarray(detection,float)
    if x.shape!=det.shape:raise ValueError('Mismatched expression/detection')
    count=np.isfinite(x).sum(axis=0);eligible=count>=3
    mean=np.full(x.shape[1],np.nan);mdet=np.full(x.shape[1],np.nan)
    np.divide(np.nansum(x,axis=0),count,out=mean,where=count>0)
    np.divide(np.nansum(det,axis=0),count,out=mdet,where=count>0)
    eligible &= np.isfinite(mean);result={'mean':mean,'mean_detection':mdet,'eligible':eligible,'top':[],'runner':[],'frequency':np.full(x.shape[1],np.nan),'valid':0,'status':'NOT_EVALUABLE','reason':'INSUFFICIENT_CATEGORIES_OR_SIGNAL'}
    if eligible.sum()<2 or np.nanmax(mdet[eligible])<.01:return result
    mx=mean[eligible].max();top=np.flatnonzero(eligible&np.isclose(mean,mx,rtol=1e-10,atol=1e-12));remaining=eligible.copy();remaining[top]=False
    runner=np.flatnonzero(remaining&np.isclose(mean,np.nanmax(mean[remaining]),rtol=1e-10,atol=1e-12)) if remaining.any() else np.array([],int)
    w=np.random.default_rng(seed).multinomial(len(x),np.full(len(x),1/len(x)),size=B)
    den=w@np.isfinite(x).astype(float);num=w@np.nan_to_num(x);boot=np.full(den.shape,np.nan);np.divide(num,den,out=boot,where=den>0)
    # Counts are bootstrap weights, not independent donor counts. Require >=3 weighted
    # eligible labels as in the original descriptive estimator; <2 types invalid.
    boot[(den<3)|~eligible[None,:]]=np.nan;ok=(np.isfinite(boot).sum(axis=1)>=2)
    freq=np.full(x.shape[1],np.nan)
    if ok.any():
        b=boot[ok];wins=np.isclose(b,np.nanmax(b,axis=1)[:,None],rtol=1e-10,atol=1e-12);share=wins/wins.sum(axis=1)[:,None];freq=share.mean(axis=0);freq[~eligible]=np.nan
    result.update(top=top.tolist(),runner=runner.tolist(),frequency=freq,valid=int(ok.sum()),status='DONE',reason='Descriptive donor-equal ranking;not differential expression')
    return result
