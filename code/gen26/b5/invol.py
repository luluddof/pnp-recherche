# Leaf S2.3.a2: search osr-r A on K^k (x) K^k with Theta = A Swap an involution and tr Theta = t > k.
# convention: A[(x,y),(u,v)] = sum_l E_l[x,u] C_l[y,v]; Swap[(x,y),(u,v)] = [x=v][y=u].
import numpy as np, sys
from scipy.optimize import least_squares
k=int(sys.argv[1]); r=int(sys.argv[2]); t=int(sys.argv[3]); trials=int(sys.argv[4]) if len(sys.argv)>4 else 20
n=k*k
S=np.zeros((n,n))
for x in range(k):
  for y in range(k):
    S[x*k+y, y*k+x]=1
def build(z):
  E=z[:r*k*k].reshape(r,k,k); C=z[r*k*k:].reshape(r,k,k)
  A=sum(np.kron(E[l],C[l]) for l in range(r))
  return A@S
def res(z):
  T=build(z)
  return np.concatenate([(T@T-np.eye(n)).ravel(), [np.trace(T)-t]])
rng=np.random.default_rng(int(sys.argv[5]) if len(sys.argv)>5 else 0)
best=1e9
for tr in range(trials):
  z0=rng.normal(size=2*r*k*k)
  sol=least_squares(res,z0,method=('lm' if n*n+1>=2*r*k*k else 'trf'),max_nfev=3000,xtol=1e-15,ftol=1e-15)
  c=np.sum(sol.fun**2)
  best=min(best,c)
  if c<1e-20:
    T=build(sol.x); ev=np.linalg.eigvals(T)
    print("FOUND trial",tr,"cost",c,"tr",np.trace(T),"maxcoef",np.abs(sol.x).max()); np.save(f"invol_{k}_{r}_{t}.npy",sol.x); break
print("k",k,"r",r,"t",t,"best cost",best)
