"""Exact absolute-tail Spearman test with weighted index permutations.

Utility only: does not authorize analysis of unadmitted source matrices.
"""
import itertools,math

def centered_ranks2(values):
    if len(values)<2 or not all(math.isfinite(x) for x in values):
        raise ValueError('At least two finite observations required')
    # Twice average ranks are integers, including ties.
    r=[sum(z<v for z in values)*2+sum(z==v for z in values)+1 for v in values]
    n=len(r)
    return [n*x-sum(r) for x in r]

def rho(x,y):
    if len(x)!=len(y):raise ValueError('Lengths differ')
    a,b=centered_ranks2(x),centered_ranks2(y)
    den=math.sqrt(sum(v*v for v in a)*sum(v*v for v in b))
    if not den:raise ValueError('Constant vector')
    return sum(u*v for u,v in zip(a,b))/den

def exact_test(x,y):
    if len(x)!=6 or len(y)!=6:raise ValueError('This locked utility requires n=6')
    a,b=centered_ranks2(x),centered_ranks2(y)
    observed=abs(sum(u*v for u,v in zip(a,b)))
    effect=rho(x,y)
    extreme=sum(abs(sum(a[i]*b[j]for i,j in enumerate(order)))>=observed
                for order in itertools.permutations(range(6)))
    leave=[];invalid=0
    for i in range(6):
        try:leave.append(rho(x[:i]+x[i+1:],y[:i]+y[i+1:]))
        except ValueError:invalid+=1
    return dict(rho=effect,p=extreme/720,extreme=extreme,permutations=720,
                loo_min=min(leave)if leave else None,loo_max=max(leave)if leave else None,
                loo_valid=len(leave),loo_invalid_constant=invalid,
                loo_sign_flip=any(v*effect<0 for v in leave))

def bh(ps):
    n=len(ps);q=[None]*n;last=1.
    for position,index in reversed(list(enumerate(sorted(range(n),key=lambda i:ps[i]),1))):
        if not 0<=ps[index]<=1:raise ValueError('Invalid P')
        last=min(last,ps[index]*n/position);q[index]=last
    return q
