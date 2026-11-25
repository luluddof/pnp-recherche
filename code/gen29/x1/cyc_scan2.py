# P_3 ring frontier with warm starts: for each bond triple descend p using the previous solution + random starts;
# also seed from embedded solutions of smaller bond triples (zero padding).
import sys, itertools, numpy as np, os
from scipy.optimize import minimize
from vp import Cyc3, fun
N=int(sys.argv[1]); bmax=int(sys.argv[2]); starts=int(os.environ.get('STARTS','6'))
rng=np.random.default_rng(1)
def embed(zsmall, bs, bb):
    ps=Cyc3(N,bs).parts(np.array(zsmall)); out=[]
    big=Cyc3(N,bb)
    for f,sh in zip(ps,big.sh):
        g=np.zeros(sh); g[:,:,:f.shape[2],:f.shape[3]]=f; out.append(g.ravel())
    return np.concatenate(out)
def opt(prob,p,z0):
    r=minimize(fun,z0,args=(prob,p),jac=True,method='L-BFGS-B',options={'maxiter':4000,'gtol':1e-14,'ftol':1e-16})
    return r.fun,r.x
sols={}  # bonds -> {p: z}
trip=sorted({min([b[i:]+b[:i] for i in range(3)]) for b in itertools.product(range(1,bmax+1),repeat=3)},key=lambda b:(np.prod(b),b))
for b in trip:
    prob=Cyc3(N,b); seeds=[]
    for bs,d in sols.items():
        if all(bs[i]<=b[i] for i in range(3)): seeds+= [(pp,embed(z,bs,b)) for pp,z in d.items()]
    d={}; p=N**3; z=None
    lo=min([pp for pp,_ in seeds],default=N**3)
    for p in range(N**3,-1,-1):
        cands=[zz+1e-3*rng.normal(size=prob.n) for pp,zz in seeds if pp<=p]
        if z is not None: cands.append(z+1e-3*rng.normal(size=prob.n))
        cands+= [rng.normal(size=prob.n)*0.5 for _ in range(starts)]
        best=(np.inf,None)
        for z0 in cands:
            f,zz=opt(prob,p,z0)
            if f<best[0]: best=(f,zz)
            if f<1e-20: break
        if best[0]>1e-9: fail=best[0]; break
        z=best[1]; d[p]=z
    sols[b]=d
    print(f"N={N} bonds={b} q={np.prod(b)} minp={min(d)} (next p tail {fail:.2e})",flush=True)
