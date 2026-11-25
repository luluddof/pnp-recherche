# independent check of Lemma 1 (compression), Lemma 2 (norm, odd traces), on optima + random MC + random BN
import numpy as np, sympy
from common import *
def analyse(N,A,Bm,Q=None,name='',exact_p=None):
    S=swap(N); K=Kmat(A,Bm); Th=K@S; n=N*N
    if Q is None: Q=sum(np.linalg.norm(a,2)*np.linalg.norm(b,2) for a,b in zip(A,Bm))
    p=np.linalg.matrix_rank(np.eye(n)-Th,tol=1e-8) if exact_p is None else exact_p
    # kernel basis of 1-Theta
    U,s,Vt=np.linalg.svd(np.eye(n)-Th); V=Vt[np.sum(s>1e-8):].T; m=V.shape[1]
    Ts=(Th+Th.T)/2
    comp_err=np.abs(V.T@Ts@V-np.eye(m)).max() if m else 0
    ev=np.linalg.eigvalsh(Ts); nge=int(np.sum(ev>=1-1e-8))
    tr=[np.trace(np.linalg.matrix_power(Ts,j)) for j in (1,3,5,7,9)]
    trok=all(abs(t)<=Q**j*N+1e-6 for t,j in zip(tr,(1,3,5,7,9)))
    return dict(name=name,N=N,q=len(A),Q=round(Q,3),p=p,m=n-p,m_svd=m,comp_err=comp_err,nge1=nge,lam_ok=nge>=n-p,
                normTs=round(ev.max() if abs(ev.max())>abs(ev.min()) else ev.min(),3),norm_ok=max(abs(ev))<=Q+1e-9,trok=trok,
                maxratio=max(abs(t)/(Q**j*N) for t,j in zip(tr,(1,3,5,7,9))))
if __name__=='__main__':
    for f,pp in OPT:
        N,A,Bm=load(f)
        pex=sympy.Matrix((swap(N)-Kmat(A,Bm)).astype(int).tolist()).rank()
        print('OPT',analyse(N,A,Bm,name=f.split('/')[-1],exact_p=pex),'stored p',pp)
    rng=np.random.default_rng(1); bad=0; cnt=0
    for t in range(400):
        N=int(rng.integers(2,7)); q=int(rng.integers(1,6))
        A=[rand_mono(N,rng,dens=rng.random()) for _ in range(q)]; Bm=[rand_mono(N,rng,dens=rng.random()) for _ in range(q)]
        r=analyse(N,A,Bm); cnt+=1
        if not(r['lam_ok'] and r['norm_ok'] and r['trok'] and r['comp_err']<1e-7): bad+=1; print('BAD MC',r)
    print('random MC',cnt,'bad',bad)
    # random BN: real gaussian factors, also near-structured (G' = G^{-T} etc)
    bad=0
    for t in range(300):
        N=int(rng.integers(2,6)); q=int(rng.integers(1,4))
        A=[rng.standard_normal((N,N)) for _ in range(q)]
        Bm=[np.linalg.inv(a).T if rng.random()<0.5 else rng.standard_normal((N,N)) for a in A]
        r=analyse(N,A,Bm)
        if not(r['lam_ok'] and r['norm_ok'] and r['trok'] and r['comp_err']<1e-6): bad+=1; print('BAD BN',r)
    print('random BN bad',bad)
