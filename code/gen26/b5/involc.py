# complex version of invol.py
import numpy as np, sys
from scipy.optimize import least_squares
k=int(sys.argv[1]); r=int(sys.argv[2]); t=int(sys.argv[3]); trials=int(sys.argv[4]); seed=int(sys.argv[5])
n=k*k; m=r*k*k
S=np.zeros((n,n))
for x in range(k):
  for y in range(k): S[x*k+y,y*k+x]=1
def build(z):
  w=z[:2*m]+1j*z[2*m:]
  E=w[:m].reshape(r,k,k); C=w[m:].reshape(r,k,k)
  return sum(np.kron(E[l],C[l]) for l in range(r))@S
def res(z):
  T=build(z); v=np.concatenate([(T@T-np.eye(n)).ravel(),[np.trace(T)-t]])
  return np.concatenate([v.real,v.imag])
rng=np.random.default_rng(seed); best=1e9
for tr in range(trials):
  sol=least_squares(res,rng.normal(size=4*m),method='trf',max_nfev=3000,xtol=1e-15,ftol=1e-15,gtol=1e-15)
  c=np.sum(sol.fun**2); best=min(best,c)
  print(tr,c,np.abs(sol.x).max(),flush=True)
  if c<1e-18: print("FOUND"); np.save(f"involc_{k}_{r}_{t}.npy",sol.x); break
print("k",k,"r",r,"t",t,"best",best)
