# progressive rounding of a hada.py solution (N, R=1): fix entries to nearest multiple of 1/den (incl. 0), re-optimise the rest.
import numpy as np, sys, time
from scipy.optimize import minimize
import hada as h
fn=sys.argv[1]; N=int(sys.argv[2]); p=int(sys.argv[3]); den=int(sys.argv[4]); batch=int(sys.argv[5]) if len(sys.argv)>5 else 8
R=1; perm=h.setup(N); z=np.load(fn); n=z.size
fixed=np.zeros(n,bool)
def opt(z,fixed,its=3000):
    free=~fixed
    def f(zf):
        zz=z.copy(); zz[free]=zf; v,g=h.fg(zz,N,R,p,perm); return v,g[free]
    r=minimize(f,z[free],jac=True,method='L-BFGS-B',options=dict(maxiter=its,maxcor=50,gtol=1e-15,ftol=1e-22))
    zz=z.copy(); zz[free]=r.x; return zz,r.fun
z,f=opt(z,fixed); print('start',f,flush=True)
while (~fixed).any():
    free=np.where(~fixed)[0]; tgt=np.round(z[free]*den)/den; d=np.abs(z[free]-tgt)
    o=np.argsort(d)[:batch]; zt=z.copy(); zt[free[o]]=tgt[o]; fx=fixed.copy(); fx[free[o]]=True
    z2,f2=opt(zt,fx)
    if f2<1e-16: z,fixed=z2,fx; print(f"fixed {fixed.sum()}/{n} tail {f2:.1e} maxdev {d[o].max():.3f}",flush=True)
    else:
        if batch>1: batch=max(1,batch//2); print('batch->',batch,'tail',f2,flush=True); continue
        print('STUCK: cannot round entry, tail',f2,'dev',d[o].max(),flush=True); break
np.save(fn.replace('.npy',f'_round{den}.npy'),z); print('done fixed',fixed.sum(),'of',n)
print('rank sum',sum(h.ranks(z,N,R,perm)))
