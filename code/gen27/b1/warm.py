import numpy as np, sys, hada as h
from scipy.optimize import minimize
fn,N,p=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]); perm=h.setup(N); z=np.load(fn)
r=minimize(h.fg,z,args=(N,1,p,perm),jac=True,method='L-BFGS-B',options=dict(maxiter=8000,maxcor=50,gtol=1e-16,ftol=1e-22))
X=r.x[:N*N].reshape(1,N,N); Y=r.x[N*N:].reshape(1,N,N)
allsv=np.sort(np.concatenate([np.linalg.svd(M,compute_uv=False) for M in h.blocks(X,Y,perm)]))[::-1]
print(fn,'target',p,'tail',r.fun,'sv[p-1],sv[p],sv[p+1]',allsv[p-1],allsv[p],allsv[p+1],'nit',r.nit,'maxabs',np.abs(r.x).max())
np.save(fn.replace('.npy',f'_w{p}.npy'),r.x)
