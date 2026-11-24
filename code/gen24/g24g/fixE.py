# [N] For FIXED E_1..E_r, minimise over C_l the tail sum_{i>p} sigma_i(Swap - sum E_l (x) C_l)^2 (L-BFGS, complex C).
# Tests whether a non-triangularizable E-span (e.g. span(I, D, w0)) admits C's with rank below N(N-1)/2.
import numpy as np, sys, time
from scipy.optimize import minimize
def swap(N):
    S=np.zeros((N*N,N*N))
    for x in range(N):
        for y in range(N): S[x*N+y, y*N+x]=1
    return S
def run(N, p, Es, R=4, seed=0, maxiter=4000):
    r=len(Es); S=swap(N).astype(complex); rng=np.random.default_rng(seed); best=1e9
    def fg(z):
        w=z[:len(z)//2]+1j*z[len(z)//2:]; C=w.reshape(r,N,N)
        K=sum(np.kron(E,c) for E,c in zip(Es,C)); M=S-K
        U,s,Vh=np.linalg.svd(M); Et=(U[:,p:]*s[p:])@Vh[p:,:]
        f=float(np.sum(s[p:]**2)); Et=Et.reshape(N,N,N,N)   # (x,y,u,v)
        g=np.stack([-2*np.einsum('xyuv,xu->yv',Et,E.conj()) for E in Es])
        g=g.ravel(); return f, np.concatenate([g.real,g.imag])
    for t in range(R):
        z0=rng.standard_normal(2*r*N*N)*rng.choice([.3,1.,2.])
        res=minimize(fg,z0,jac=True,method='L-BFGS-B',options=dict(maxiter=maxiter,ftol=1e-16,gtol=1e-12,maxcor=30))
        best=min(best,res.fun)
    return best
if __name__=="__main__":
    N=int(sys.argv[1]); kind=sys.argv[2]; ps=list(map(int,sys.argv[3].split(',')))
    rng=np.random.default_rng(7)
    I=np.eye(N); D=np.diag(np.arange(1,N+1)*1.0); w0=np.fliplr(np.eye(N))
    if kind=='IDw0': Es=[I,D,w0]
    elif kind=='rand3': Es=[rng.standard_normal((N,N)) for _ in range(3)]
    elif kind=='tri3': Es=[np.triu(rng.standard_normal((N,N))) for _ in range(3)]
    elif kind=='Iw0': Es=[I,w0]
    for p in ps:
        t0=time.time(); b=run(N,p,Es)
        print(f"N={N} E={kind} p={p} best_tail={b:.6f} ({time.time()-t0:.0f}s)",flush=True)
