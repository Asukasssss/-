"""Pure statistical kernels for versioned COAD patient analysis."""
import hashlib
import numpy as np
from scipy import stats


def stable_seed(text):return int(hashlib.sha256(text.encode()).hexdigest()[:8],16)


def bh(p):
    p=np.asarray(p,dtype=float)
    assert np.isfinite(p).all() and ((p>=0)&(p<=1)).all()
    order=np.argsort(p);n=len(p)
    adjusted=np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1]
    out=np.empty(n);out[order]=np.minimum(adjusted,1)
    return out


def row_rho(x,y):
    x=stats.rankdata(x,axis=-1,method='average');y=stats.rankdata(y,axis=-1,method='average')
    x=x-x.mean(axis=-1,keepdims=True);y=y-y.mean(axis=-1,keepdims=True)
    with np.errstate(invalid='ignore',divide='ignore'):
        return np.sum(x*y,axis=-1)/np.sqrt(np.sum(x*x,axis=-1)*np.sum(y*y,axis=-1))


def association(x,y,seed,permutations=9999,bootstrap=4000,minimum_n=8,loo=False):
    good=np.isfinite(x)&np.isfinite(y);x=np.asarray(x)[good];y=np.asarray(y)[good]
    n=len(x);base=dict(n=n,n_unique_metabolite=len(np.unique(x)),n_unique_rna=len(np.unique(y)))
    if n<minimum_n:return dict(base,status='NOT_EVALUABLE',reason='N_BELOW_MINIMUM')
    if np.ptp(x)==0 or np.ptp(y)==0:return dict(base,status='NOT_EVALUABLE',reason='CONSTANT_VALUES')
    rho=float(row_rho(x,y));assert np.isclose(rho,stats.spearmanr(x,y).statistic,atol=1e-12)
    rx=stats.rankdata(x);ry=stats.rankdata(y);rx-=rx.mean();ry-=ry.mean()
    rng=np.random.default_rng(seed)
    perm=np.array([rng.permutation(n) for _ in range(permutations)])
    null=np.einsum('ij,j->i',ry[perm],rx)/np.sqrt(np.dot(rx,rx)*np.dot(ry,ry))
    p=(int(np.count_nonzero(np.abs(null)>=abs(rho)-1e-12))+1)/(permutations+1)
    draws=rng.integers(0,n,size=(bootstrap,n));boot=row_rho(x[draws],y[draws]);boot=boot[np.isfinite(boot)]
    ci=np.quantile(boot,[.025,.975]) if len(boot) else [None,None]
    result=dict(base,status='DONE',reason='COMPUTED',effect=rho,p_value=p,
       ci_lower=float(ci[0]) if ci[0] is not None else None,ci_upper=float(ci[1]) if ci[1] is not None else None,
       bootstrap_valid=len(boot),permutations=permutations,random_seed=seed)
    if loo:
        ix=np.broadcast_to(np.arange(n),(n,n));ix=ix[~np.eye(n,dtype=bool)].reshape(n,n-1)
        rr=row_rho(x[ix],y[ix]);rr=rr[np.isfinite(rr)]
        result.update(loo_valid=len(rr),loo_min=float(rr.min()) if len(rr) else None,
          loo_max=float(rr.max()) if len(rr) else None,loo_max_abs_delta=float(np.max(np.abs(rr-rho))) if len(rr) else None,
          loo_any_sign_change=bool((rr*rho<0).any()))
    return result


def paired_rna(tumor,normal,seed,bootstrap=4000,minimum_n=8):
    good=np.isfinite(tumor)&np.isfinite(normal);d=np.asarray(tumor)[good]-np.asarray(normal)[good]
    base=dict(n=len(d),n_reference=len(d),n_nonzero_differences=int(np.count_nonzero(d)))
    if len(d)<minimum_n:return dict(base,status='NOT_EVALUABLE',reason='N_BELOW_MINIMUM')
    p=1.0 if not np.any(d) else float(stats.wilcoxon(d,zero_method='wilcox',correction=False,alternative='two-sided',method='approx').pvalue)
    rng=np.random.default_rng(seed);draws=rng.integers(0,len(d),size=(bootstrap,len(d)))
    boot=np.median(d[draws],axis=1);lo,hi=np.quantile(boot,[.025,.975])
    return dict(base,status='DONE',reason='COMPUTED',effect=float(np.median(d)),p_value=p,
      ci_lower=float(lo),ci_upper=float(hi),bootstrap_valid=len(boot),random_seed=seed)
