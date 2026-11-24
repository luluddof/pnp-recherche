import numpy as np, sys
from scipy.optimize import least_squares
import crossing_ap as ca
# re-implement with access to parameters: long runs for N=4, (r,c)=(3,3)
N=int(sys.argv[1]); r=int(sys.argv[2]); c=int(sys.argv[3]); seed=int(sys.argv[4]); nfev=int(sys.argv[5])
n=N*N; rng=np.random.default_rng(seed)
T=np.zeros((N,)*4)
for x in range(N):
    for y in range(N): T[x,y,y,x]=1
Trow=T.reshape(n,n); perm=np.arange(n*n).reshape(N,N,N,N).transpose(0,2,1,3).reshape(-1); In=np.eye(n)
def unpack(p):
    i=0; out=[]
    for k in (r,r,c,c): out.append(p[i:i+n*k].reshape(n,k)); i+=n*k
    return out
def res(p):
    U,V,X,Y=unpack(p); Rm=(U@V.T).ravel()-Trow.ravel(); Rm[perm]+=(X@Y.T).ravel(); return Rm
def jac(p):
    U,V,X,Y=unpack(p)
    JX=np.zeros((n*n,n*c)); JY=np.zeros((n*n,n*c)); JX[perm]=np.kron(In,Y); JY[perm]=np.einsum('ik,jl->ijlk',X,In).reshape(n*n,n*c)
    return np.hstack([np.kron(In,V), np.einsum('ik,jl->ijlk',U,In).reshape(n*n,n*r), JX, JY])
for t in range(6):
    p0=rng.normal(size=2*n*(r+c))*0.5
    s=least_squares(res,p0,jac=jac,method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=nfev)
    U,V,X,Y=unpack(s.x); A=U@V.T; B=X@Y.T
    print(f"try {t}: res^2={np.sum(s.fun**2):.3e} nfev={s.nfev} |A|_F={np.linalg.norm(A):.2f} |B|_F={np.linalg.norm(B):.2f} sv(A)={np.round(np.linalg.svd(A,compute_uv=False)[:r],3)} sv(B)={np.round(np.linalg.svd(B,compute_uv=False)[:c],3)}",flush=True)
    if np.sum(s.fun**2)<1e-3: np.save(f"sol_N{N}_{r}_{c}_s{seed}_t{t}.npy", s.x)
