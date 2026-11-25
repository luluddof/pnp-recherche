import sys, math; sys.path.insert(0,'.')
from model2 import *
from scipy.optimize import minimize_scalar, minimize
TAUS=((1.0,'2^n'),(0.386,'PPSZ'))
def best_eps(lam, **kw):
    # optimise eps for the paper-type schedule (u_d = eps lam/2^d)
    f=lambda le: run(lam, paper_sched(lam, math.exp(le)), **kw)
    res=minimize_scalar(f, bounds=(math.log(1e-6), math.log(0.2)), method='bounded', options={'xatol':1e-5})
    # coarse scan to avoid local minima from floors/infeasibility
    bestv=res.fun; beste=math.exp(res.x)
    for k in range(200):
        e=10**(-6+k*0.025)
        v=f(math.log(e))
        if v<bestv: bestv, beste=v, e
    return bestv, beste
rows=[
 ("A  paper as written", lambda L: run(L, paper_sched(L,2**-17))),
 ("B  + eps at legal max of (D.6)", lambda L: run(L, paper_sched(L,eps_legal(L)[0]))),
 ("C1 + Lemma D.5 budget = its proved loss (8hD->6hD, 8M(a+1)->q(a+1))", lambda L: run(L, paper_sched(L,eps_legal(L)[0]),kappa='18')),
 ("C2 + minimal multiplicity step h=ceil(M/(u-1)) (='LP taper')", lambda L: run(L, paper_sched(L,eps_legal(L)[0]),kappa='lin')),
 ("D  + packing slack (nu=2P(1-1/r), child 2PM/r), eps re-optimised", lambda L: best_eps(L,kappa='lin',pack='slack')[0]),
 ("E  + incidence network M0=km (w 45.1->12.78)", lambda L: best_eps(L,kappa='lin',pack='slack',w=W_INC)[0]),
 ("F  + adaptive degree bound Delta_d=max(D0,2L_parent)", lambda L: best_eps(L,kappa='lin',pack='slack',w=W_INC,Dmode='adapt')[0]),
]
res={}
for nm,fn in rows:
    out=[]
    for tau,tn in TAUS:
        th=threshold(lambda L: fn(L), tau)
        out.append(th)
    res[nm]=out
    print(f"{nm:75s}  vs 2^n: {out[0]:.3g}   vs PPSZ: {out[1]:.3g}   E*sqrt(lam)/n @1e8: {fn(1e8)*1e4:.1f}", flush=True)
