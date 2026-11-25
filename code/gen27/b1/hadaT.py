# hada with tied Y = c * X^T (c fixed 1/2 after scaling X)  -> symmetric-type sub-ansatz
import numpy as np, sys, time, hada as h
from scipy.optimize import minimize
N=int(sys.argv[1]); ps=list(map(int,sys.argv[2].split(','))); seeds=int(sys.argv[3]); perm=h.setup(N)
def fg(x,p):
    X=x.reshape(N,N); z=np.concatenate([X.ravel(),0.5*X.T.ravel()])
    f,g=h.fg(z,N,1,p,perm); gX=g[:N*N].reshape(N,N)+0.5*g[N*N:].reshape(N,N).T
    return f,gX.ravel()
for p in ps:
    for sd in range(seeds):
        x=np.random.default_rng(sd).standard_normal(N*N); t=time.time()
        r=minimize(fg,x,args=(p,),jac=True,method='L-BFGS-B',options=dict(maxiter=4000,maxcor=50,gtol=1e-15,ftol=1e-20))
        z=np.concatenate([r.x,0.5*r.x.reshape(N,N).T.ravel()]); rk=h.ranks(z,N,1,perm)
        print(f"tied N={N} p={p} seed={sd}: tail={r.fun:.2e} ranksum={sum(rk)} {rk} t={time.time()-t:.0f}s",flush=True)
        if r.fun<1e-14: np.save(f"hT_N{N}_p{p}_s{sd}.npy",z)
