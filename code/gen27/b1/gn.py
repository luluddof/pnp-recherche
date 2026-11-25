# Gauss-Newton polish: residual = M_kappa V_kappa (V_kappa = current right null-space estimate, sizes from global tail split)
import numpy as np, sys, hada as h
from scipy.optimize import least_squares
fn,N,p=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]); perm=h.setup(N); z=np.load(fn); n=N*N
def split(z):
    X=z[:n].reshape(1,N,N); Y=z[n:].reshape(1,N,N); svd=[np.linalg.svd(M) for M in h.blocks(X,Y,perm)]
    allS=np.concatenate([s[1] for s in svd]); o=np.argsort(allS); mask=np.zeros(n,bool); mask[o[:n-p]]=True
    Vs=[svd[b][2][mask[b*N:(b+1)*N]].T for b in range(N)]; return Vs, np.sort(allS)[::-1]
for it in range(12):
    Vs,allS=split(z)
    def res(zz):
        X=zz[:n].reshape(1,N,N); Y=zz[n:].reshape(1,N,N)
        return np.concatenate([(M@V).ravel() for M,V in zip(h.blocks(X,Y,perm),Vs)])
    r=least_squares(res,z,method='lm',max_nfev=200*len(z)//len(z)*20,xtol=1e-15,ftol=1e-15,gtol=1e-15)
    z=r.x; _,allS=split(z)
    print(it,'cost',r.cost,'sv[p-1],sv[p],sv[p+1]',allS[p-1],allS[p],allS[p+1],'maxabs',np.abs(z).max(),flush=True)
np.save(fn.replace('.npy',f'_gn{p}.npy'),z)
