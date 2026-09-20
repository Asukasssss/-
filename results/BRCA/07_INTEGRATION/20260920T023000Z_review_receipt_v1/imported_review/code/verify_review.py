"""Review-artifact integrity and synthetic algebra checks; not a patient reanalysis."""
from pathlib import Path
import csv, json
import numpy as np
from scipy import stats

def read(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def main():
    root=Path(__file__).resolve().parents[1]
    a=read(root/"tables/association_excerpts_14.tsv")
    s=read(root/"tables/source_stability_excerpts_17.tsv")
    assert len(a)==14 and len({r["relation_id"] for r in a})==14
    assert len(s)==17 and len({r["gene"] for r in s})==17
    for r in a:
        for stem in ("ER_PC12","ER_MonoEndoFib"):
            rho=float(r[stem+"_rho"])
            assert -1 <= rho <= 1
            assert float(r[stem+"_ci_lower"]) <= rho <= float(r[stem+"_ci_upper"])
            for suffix in ("_p","_q"):
                assert 0 <= float(r[stem+suffix]) <= 1
    for r in s:
        stable=(r["Wu_top"]==r["Pal_top"] and min(float(r["Wu_bootstrap_top_frequency"]),
                 float(r["Pal_bootstrap_top_frequency"]))>=.8)
        assert stable==(r["two_source_stable_rule"]=="True")
    rng=np.random.default_rng(20260919)
    e1=e2=0.
    for _ in range(30):
        er=np.r_[np.zeros(30),np.ones(31)]
        cov=rng.normal(size=(61,3))
        v=np.column_stack([er]+[stats.rankdata(cov[:,j]) for j in range(3)])
        v=(v-v.mean(0))/v.std(0)
        z=np.column_stack([np.ones(61),v]); q=np.linalg.qr(z,mode="reduced")[0]
        x=stats.rankdata(rng.normal(size=61)); y=stats.rankdata(rng.normal(size=61))
        ex=x-q@(q.T@x); ey=y-q@(q.T@y)
        rho=ex@ey/(np.linalg.norm(ex)*np.linalg.norm(ey))
        ax=x-z@np.linalg.lstsq(z,x,rcond=None)[0]
        ay=y-z@np.linalg.lstsq(z,y,rcond=None)[0]
        e1=max(e1,abs(rho-np.corrcoef(ax,ay)[0,1]))
        for _ in range(20):
            per=ey.copy()
            for level in (0,1):
                ix=np.flatnonzero(er==level)
                per[ix]=rng.permutation(ey[ix])
            fast=per@ex/(np.linalg.norm(ex)*np.sqrt(np.sum(per**2)-np.sum((per@q)**2)))
            ystar=q@(q.T@y)+per
            yp=ystar-q@(q.T@ystar)
            slow=ex@yp/(np.linalg.norm(ex)*np.linalg.norm(yp))
            e2=max(e2,abs(fast-slow))
    assert e1<1e-12 and e2<1e-12
    print(json.dumps({"excerpt_checks":"PASS","synthetic_QR_OLS_max_error":e1,
       "synthetic_permutation_identity_max_error":e2,"patient_tests":0,
       "scope":"Not full348 BH,not full117 biological validation"},ensure_ascii=False,indent=2))
if __name__=="__main__":
    main()
