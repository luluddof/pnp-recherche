# Greedy (OMP-like) support search for the Bell model, scored by rw-LAD; exact verification of frontier points.
import numpy as np, sys, json, time
from bell import *
k=int(sys.argv[1]); qmax=int(sys.argv[2]); add=int(sys.argv[3]) if len(sys.argv)>3 else 1; seed=int(sys.argv[4]) if len(sys.argv)>4 else 0
rng=np.random.default_rng(seed)
b=beta(k); n=4**k
allz=np.arange(n)
S=[0]
front={}; keep={}
t0=time.time()
while len(S)<=qmax:
    M=charmat(k,S); e,x=rw_lad(M,b,iters=4)
    p,nz=exact_errors(M,b,x)
    if nz not in front or p<front[nz]: front[nz]=p; keep[nz]=[(int(s),round(float(v),4)) for s,v in zip(S,x) if abs(v)>1e-9]
    r=M@x-b; E=np.abs(r)>1e-7
    sig=np.where(E,np.sign(r),0.0)
    sc=np.abs(fwht(sig)); sc[S]=-1
    sc=sc+1e-9*rng.random(n)
    new=list(np.argsort(-sc)[:add])
    S=S+[int(t) for t in new]
    print(len(S)-add, 'p_exact',p,'nnz',nz, f'{time.time()-t0:.0f}s', flush=True)
best=sorted(front.items())
print('FRONTIER k=%d'%k, best)
print('S6',keep.get(6))
