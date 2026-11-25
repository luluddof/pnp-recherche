# real variable projection for Swap = R + C, rank R <= p, C = sum_l G_l (x) H_l.  Objective: tail = sum_{i>p} sigma_i(S-K)^2
# (mode 'f2') or log-tail sum log(sigma_i^2+eps) with eps continuation (mode 'log').  Diagnostics: spectral profile of Ts.
import numpy as np, sys, time, json, argparse
from scipy.optimize import minimize
def swap(N):
    X=np.zeros((N,)*4)
    for x in range(N):
        for y in range(N): X[x,y,y,x]=1
    return X
def unpack(z,N,q): s=N*N; return z[:q*s].reshape(q,N,N), z[q*s:].reshape(q,N,N)
def fg(z,F,N,p,q,mode,eps,lam):
    s=N*N; G,H=unpack(z,N,q)
    B=np.einsum('kxu,kyv->xyuv',G,H); M=(F-B).reshape(s,s)
    U,S,Vh=np.linalg.svd(M)
    if mode=='f2': w=2*S[p:]; f=float(np.sum(S[p:]**2))
    else: w=2*S[p:]/(S[p:]**2+eps); f=float(np.sum(np.log(S[p:]**2+eps)-np.log(eps)))
    E=((U[:,p:]*w)@Vh[p:,:]).reshape(N,N,N,N)
    gG=-np.einsum('xyuv,kyv->kxu',E,H); gH=-np.einsum('xyuv,kxu->kyv',E,G)
    if lam>0:  # balance penalty ||G||^2-||H||^2 per term (gauge, does not change K)
        a=np.sum(G**2,axis=(1,2)); b=np.sum(H**2,axis=(1,2)); d=a-b
        f+=lam*float(np.sum(d**2)); gG+=lam*4*d[:,None,None]*G; gH-=lam*4*d[:,None,None]*H
    return f,np.concatenate([gG.ravel(),gH.ravel()])
def diag(z,N,p,q):
    s=N*N; G,H=unpack(z,N,q); F=swap(N)
    K=np.einsum('kxu,kyv->xyuv',G,H).reshape(s,s); Sm=F.reshape(s,s)
    sv=np.linalg.svd(Sm-K,compute_uv=False)
    Th=K@Sm; Ts=(Th+Th.T)/2; ev=np.linalg.eigvalsh(Ts)
    Q=sum(np.linalg.norm(G[l],2)*np.linalg.norm(H[l],2) for l in range(q))
    m=s-p
    out=dict(tail=float(np.sum(sv[p:]**2)),sp=float(sv[p-1]) if p>0 else 0.,sp1=float(sv[p]),gap=float(sv[p-1]/max(sv[p],1e-300)) if p>0 else 0.,
        coef=float(max(abs(G).max(),abs(H).max())),Q=float(Q),normTs=float(abs(ev).max()),trTh=float(np.trace(Th)),
        nge1=int(np.sum(ev>=1-1e-6)),m=m,nneg=int(np.sum(ev<-1e-6)),evmin=float(ev.min()),evmax=float(ev.max()),
        thev=sorted(np.linalg.eigvals(Th),key=lambda c:-abs(c))[:6], ev_sorted_low=[round(float(x),3) for x in ev[:6]],ev_sorted_high=[round(float(x),3) for x in ev[-6:]])
    out['thev']=[str(np.round(c,2)) for c in out['thev']]; return out
if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('N',type=int); ap.add_argument('p',type=int); ap.add_argument('q',type=int)
    ap.add_argument('--R',type=int,default=5); ap.add_argument('--seed',type=int,default=0); ap.add_argument('--it',type=int,default=3000)
    ap.add_argument('--mode',default='f2'); ap.add_argument('--lam',type=float,default=0.0); ap.add_argument('--warm',default=None)
    ap.add_argument('--tmax',type=float,default=500); ap.add_argument('--tag',default='')
    a=ap.parse_args(); N,p,q=a.N,a.p,a.q; F=swap(N); rng=np.random.default_rng(a.seed); best=None; t0=time.time()
    for r in range(a.R):
        if time.time()-t0>a.tmax: break
        if a.warm: z0=np.load(a.warm); z0=z0+rng.standard_normal(z0.size)*0.05*(r>0)
        else: z0=rng.standard_normal(2*q*N*N)*rng.choice([.4,.7,1.])
        if a.mode=='log':
            z=z0
            for eps in [1e-1,1e-2,1e-3,1e-4,1e-6]:
                z=minimize(lambda z: fg(z,F,N,p,q,'log',eps,a.lam),z,jac=True,method='L-BFGS-B',options=dict(maxiter=a.it//5,ftol=1e-15,gtol=1e-12,maxcor=30)).x
            z=minimize(lambda z: fg(z,F,N,p,q,'f2',0,a.lam),z,jac=True,method='L-BFGS-B',options=dict(maxiter=a.it,ftol=1e-16,gtol=1e-13,maxcor=30)).x
        else:
            z=minimize(lambda z: fg(z,F,N,p,q,'f2',0,a.lam),z0,jac=True,method='L-BFGS-B',options=dict(maxiter=a.it,ftol=1e-16,gtol=1e-13,maxcor=30)).x
        d=diag(z,N,p,q)
        print(f"r{r} tail={d['tail']:.2e} gap={d['gap']:.1e} coef={d['coef']:.2f} Q={d['Q']:.1f} ||Ts||={d['normTs']:.2f} trTh={d['trTh']:.1f} #ev>=1={d['nge1']}(m={d['m']}) #neg={d['nneg']} evmin={d['evmin']:.2f} evmax={d['evmax']:.2f} t={time.time()-t0:.0f}s",flush=True)
        if best is None or d['tail']<best[0]: best=(d['tail'],z,d)
    np.save(f"z_{N}_{p}_{q}{a.tag}.npy",best[1]); json.dump(best[2],open(f"d_{N}_{p}_{q}{a.tag}.json","w"))
    print("BEST",json.dumps(best[2]))
