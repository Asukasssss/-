"""Partial ranks and approximate Freedman-Lane with nuisance refitting."""
import numpy as np
from scipy.stats import rankdata
from patient_statistics import average_ranks

def partial(x,y,z):
    rx=rankdata(x);ry=rankdata(y)
    inv=np.linalg.pinv(z);ex=rx-z@(inv@rx);ey=ry-z@(inv@ry)
    denom=np.linalg.norm(ex)*np.linalg.norm(ey)
    if denom<1e-10:return None
    return float(np.clip(ex@ey/denom,-1,1))

def fl_null(rx,ry,z,perm):
    inv=np.linalg.pinv(z);ex=rx-z@(inv@rx);ey=ry-z@(inv@ry)
    # Adding the fitted nuisance component then projecting it out is algebraically zero.
    ep=ey[perm];ep=ep-(ep@inv.T)@z.T
    den=np.linalg.norm(ep,axis=1)*np.linalg.norm(ex)
    return (ep@ex)/den

def bootstrap_plan(z,seed,b):
    rng=np.random.default_rng(seed);n=len(z);idx=rng.integers(0,n,size=(b,n))
    zb=z[idx];ranks=np.linalg.matrix_rank(zb);valid=ranks==z.shape[1]
    return idx,zb,np.linalg.pinv(zb),valid

def adjusted(x,y,z,seed,spec,bootplan):
    n=len(x);k=z.shape[1];rho=partial(x,y,z)
    if rho is None:return {'status':'NOT_EVALUABLE','reason':'CONSTANT_RESIDUAL'}
    rx=rankdata(x);ry=rankdata(y)
    # Independent OLS implementation, not the projection used to calculate rho.
    ex=rx-z@np.linalg.lstsq(z,rx,rcond=None)[0];ey=ry-z@np.linalg.lstsq(z,ry,rcond=None)[0]
    check=float(np.corrcoef(ex,ey)[0,1]);assert abs(check-rho)<1e-10
    rng=np.random.default_rng(seed);perm=np.array([rng.permutation(n) for _ in range(spec['permutations'])])
    null=fl_null(rx,ry,z,perm);assert np.isfinite(null).all()
    p=(np.count_nonzero(np.abs(null)>=abs(rho)-1e-12)+1)/(len(null)+1)
    ix,zb,pinv,valid=bootplan
    bx=average_ranks(x[ix]);by=average_ranks(y[ix])
    bx-=np.einsum('bij,bj->bi',zb,np.einsum('bij,bj->bi',pinv,bx))
    by-=np.einsum('bij,bj->bi',zb,np.einsum('bij,bj->bi',pinv,by))
    den=np.linalg.norm(bx,axis=1)*np.linalg.norm(by,axis=1)
    use=valid & (den>1e-10)
    rr=np.sum(bx[use]*by[use],axis=1)/den[use]
    ci=np.quantile(rr,[.025,.975]) if len(rr)>=.8*len(ix) else [None,None]
    loo=[]
    for i in range(n):
        mask=np.arange(n)!=i;zz=z[mask]
        if np.linalg.matrix_rank(zz)==k and (n-1-k-1)>=spec['minimum_full_model_residual_df']:
            r=partial(x[mask],y[mask],zz)
            if r is not None:loo.append(r)
    return {'status':'DONE','reason':'APPROXIMATE_FREEDMAN_LANE','effect':rho,'p_value':float(p),
        'ci_lower':None if ci[0] is None else float(ci[0]),'ci_upper':None if ci[1] is None else float(ci[1]),
        'ci_status':'DONE' if ci[0] is not None else 'INSUFFICIENT_VALID_BOOTSTRAPS',
        'bootstrap_valid':int(use.sum()),'bootstrap_invalid':int((~use).sum()),'permutations':len(null),'random_seed':seed,
        'ols_max_abs_error':abs(check-rho),'loo_valid':len(loo),'loo_invalid':n-len(loo),
        'loo_min':min(loo) if loo else None,'loo_max':max(loo) if loo else None,
        'loo_max_abs_delta':max(abs(r-rho) for r in loo) if loo else None,
        'loo_any_sign_change':any(r*rho<0 for r in loo)}
