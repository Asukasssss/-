"""Synthetic checks only; independent direct refits validate nuisance permutation."""
import json
import numpy as np
from scipy import stats
from covariate_statistics import partial,fl_null,bootstrap_plan,adjusted
from patient_statistics import bh

def validate():
 rng=np.random.default_rng(85923);n=33
 age=rng.normal(size=n);g=np.arange(n)%2;z=np.column_stack([np.ones(n),age,g])
 x=np.round(2*age+rng.normal(size=n),1);y=1.4*age+.6*x+rng.normal(size=n)
 rx=stats.rankdata(x);ry=stats.rankdata(y);r=partial(x,y,z)
 design=np.column_stack([z,rx]);b=np.linalg.lstsq(design,ry,rcond=None)[0];e=ry-design@b;df=n-design.shape[1]
 se=np.sqrt((e@e/df)*np.linalg.inv(design.T@design)[-1,-1]);t=b[-1]/se
 assert abs(t-r*np.sqrt(df/(1-r*r)))<1e-10
 perms=np.array([rng.permutation(n) for _ in range(199)])
 null=fl_null(rx,ry,z,perms);fit=z@np.linalg.lstsq(z,ry,rcond=None)[0];res=ry-fit
 direct=[]
 for ix in perms:
  yp=fit+res[ix];ey=yp-z@np.linalg.lstsq(z,yp,rcond=None)[0];ex=rx-z@np.linalg.lstsq(z,rx,rcond=None)[0]
  direct.append(np.corrcoef(ex,ey)[0,1])
 assert np.max(np.abs(null-direct))<1e-10
 assert partial(np.ones(n),y,z) is None
 assert np.linalg.matrix_rank(np.column_stack([z,z[:,-1]]))<4
 bplan=bootstrap_plan(z,93,400)
 s={'permutations':499,'minimum_full_model_residual_df':8}
 a=adjusted(x,y,z,54,s,bplan);assert a['bootstrap_valid']>300
 p=np.array([.02,.8,.002,.1,.7]);order=np.argsort(p);q=np.empty(len(p));best=1
 for j in range(len(p)-1,-1,-1):best=min(best,p[order[j]]*len(p)/(j+1));q[order[j]]=best
 assert np.allclose(q,bh(p))
 # Monte Carlo smoke test under a null rank model with nuisance signal.
 rejected=0;strong=0
 for i in range(100):
  xx=stats.rankdata(rng.normal(size=n)+age)
  yy=z@np.array([0.,1.3,.4])+rng.normal(size=n)
  ex=xx-z@np.linalg.lstsq(z,xx,rcond=None)[0];ey=yy-z@np.linalg.lstsq(z,yy,rcond=None)[0]
  rr=np.corrcoef(ex,ey)[0,1];pp=np.array([rng.permutation(n) for _ in range(199)])
  pv=(np.count_nonzero(np.abs(fl_null(xx,yy,z,pp))>=abs(rr))+1)/200
  rejected+=pv<.05
  alternative=yy+2*ex;ee=alternative-z@np.linalg.lstsq(z,alternative,rcond=None)[0]
  ar=np.corrcoef(ex,ee)[0,1]
  ap=(np.count_nonzero(np.abs(fl_null(xx,alternative,z,pp))>=abs(ar))+1)/200
  strong+=ap<.05
 assert rejected<=15 and strong>=95,(rejected,strong)
 return {'status':'PASS','full_OLS_t_error':float(abs(t-r*np.sqrt(df/(1-r*r)))),'direct_refit_permutation_max_error':float(np.max(np.abs(null-direct))),'null_simulations':100,'null_p_lt_005':int(rejected),'strong_alternative_detected':int(strong),'bootstrap_valid_smoke':a['bootstrap_valid'],'checks':['tied ranks','independent OLS t and rho','199 identical-permutation full nuisance refits','constant residual','rank deficient design','paired bootstrap','independent BH'],'limitations':['100 null simulations are an implementation smoke check, not a proof of calibration for actual patients','Approximate exchangeable-error permutation; no heteroskedasticity correction']}
if __name__=='__main__':print(json.dumps(validate(),indent=2))
