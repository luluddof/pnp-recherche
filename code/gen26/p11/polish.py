import numpy as np, sys
from scipy.optimize import minimize
from vp11 import fg, swap
import gauge
f,N,p,q,out=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),sys.argv[5]; z=np.load(f); s=N*N; F=swap(N)
G,H=z[:q*s].reshape(q,N,N),z[q*s:].reshape(q,N,N)
sc=[np.sqrt(np.linalg.norm(h)/np.linalg.norm(g)) for g,h in zip(G,H)]; G=np.array([g*c for g,c in zip(G,sc)]); H=np.array([h/c for h,c in zip(H,sc)])
G2,H2=gauge.run(G,H,N,1); zb=np.concatenate([G2.ravel(),H2.ravel()])
zb=minimize(lambda z: fg(z,F,N,p,q,'f2',0,0),zb,jac=True,method='L-BFGS-B',options=dict(maxiter=6000,ftol=1e-30,gtol=1e-18,maxcor=50)).x
K=np.einsum('kxu,kyv->xyuv',zb[:q*s].reshape(q,N,N),zb[q*s:].reshape(q,N,N)).reshape(s,s); sv=np.linalg.svd(F.reshape(s,s)-K,compute_uv=False)
print(f"polished: Qmin={gauge.Qop(G2,H2):.2f} coef={abs(zb).max():.2f} sv[p-1]={sv[p-1]:.3e} sv[p]={sv[p]:.3e} gap={sv[p-1]/sv[p]:.2e}")
np.save(out,zb)
