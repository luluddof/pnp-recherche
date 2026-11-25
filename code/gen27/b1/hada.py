# Z2^k conjugation-orbit ansatz in Fourier form [verified == orbit3 with c=1]:
#   S - K  ~  direct sum over kappa in F2^k of  M_kappa = I - sum_r (X_r P_kappa) o (P_kappa Y_r),   P_kappa: i -> i xor kappa
#   p = sum_kappa rank(M_kappa), q <= R N.   Objective: global tail (n - p smallest singular values over all blocks).
import numpy as np, argparse, time
from scipy.optimize import minimize
def setup(N):
    idx=np.arange(N); return [idx^k for k in range(N)]   # perm arrays: (A P_k)[:,j] = A[:, j^k];  (P_k B)[i,:] = B[i^k,:]
def blocks(X,Y,perm):
    N=X.shape[1]; I=np.eye(N)
    return [I-np.einsum('rij,rij->ij',X[:,:,pk],Y[:,pk,:]) for pk in perm]
def fg(z,N,R,p,perm):
    X=z[:R*N*N].reshape(R,N,N); Y=z[R*N*N:].reshape(R,N,N)
    Ms=blocks(X,Y,perm); svd=[np.linalg.svd(M) for M in Ms]
    allS=np.concatenate([s[1] for s in svd]); n=N*N
    thr=np.sort(allS)[n-p-1] if p<n else -1   # n-p smallest are tail
    order=np.argsort(allS); tailmask=np.zeros(n,bool); tailmask[order[:n-p]]=True
    f=float(np.sum(allS[tailmask]**2)); gX=np.zeros_like(X); gY=np.zeros_like(Y)
    for b,(pk,(U,S,Vh)) in enumerate(zip(perm,svd)):
        m=tailmask[b*N:(b+1)*N]
        if not m.any(): continue
        Gm=(U[:,m]*(2*S[m]))@Vh[m,:]          # df/dM
        # M = I - sum_r X[:, pk] o Y[pk, :]
        gXp=-Gm[None]*Y[:,pk,:]; gYp=-Gm[None]*X[:,:,pk]
        gX[:,:,pk]+=gXp; gY[:,pk,:]+=gYp
    return f,np.concatenate([gX.ravel(),gY.ravel()])
def ranks(z,N,R,perm,tol=1e-7):
    X=z[:R*N*N].reshape(R,N,N); Y=z[R*N*N:].reshape(R,N,N)
    out=[]
    for M in blocks(X,Y,perm):
        s=np.linalg.svd(M,compute_uv=False); out.append(int(np.sum(s>tol*max(1,s[0]))))
    return out
if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--N',type=int); ap.add_argument('--R',type=int,default=1)
    ap.add_argument('--p',type=int,nargs='+'); ap.add_argument('--seeds',type=int,default=4); ap.add_argument('--its',type=int,default=5000)
    ap.add_argument('--save',default=None); ap.add_argument('--scale',type=float,default=1.0)
    a=ap.parse_args(); N,R=a.N,a.R; perm=setup(N)
    for p in a.p:
        for sd in range(a.seeds):
            rng=np.random.default_rng(sd); t=time.time()
            z=rng.standard_normal(2*R*N*N)*a.scale/np.sqrt(N)
            res=minimize(fg,z,args=(N,R,p,perm),jac=True,method='L-BFGS-B',options=dict(maxiter=a.its,maxcor=50,gtol=1e-15,ftol=1e-20))
            rk=ranks(res.x,N,R,perm)
            print(f"Z2^k N={N} R={R} (q<={R*N}) p={p} seed={sd}: tail={res.fun:.3e} block ranks(tol1e-7) sum={sum(rk)} {rk} nit={res.nit} t={time.time()-t:.0f}s",flush=True)
            if a.save and res.fun<1e-12: np.save(f"{a.save}_N{N}_R{R}_p{p}_s{sd}.npy",res.x)
