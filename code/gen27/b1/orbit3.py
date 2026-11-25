# fast plain orbit ansatz: K = sum_{g,r} c_{g,r} (L_g A_r P_g) (x) (L_g B_r P'_g)   (matmul implementation)
import numpy as np, argparse, time
from scipy.optimize import minimize
from orbit import group, perms, swap
def build(z,L,P,P2,R):
    g,n=L.shape[0],L.shape[1]
    A=z[:R*n*n].reshape(R,n,n); B=z[R*n*n:2*R*n*n].reshape(R,n,n); c=z[2*R*n*n:].reshape(g,R)
    G=np.einsum('gik,rkl,glj->grij',L,A,P,optimize=True); H=np.einsum('gik,rkl,glj->grij',L,B,P2,optimize=True)
    return A,B,c,G,H
def Kmat(c,G,H,n):  # returns K as (xy),(uv) matrix
    Gc=(G*c[:,:,None,None]).reshape(-1,n*n); Hm=H.reshape(-1,n*n)
    KR=Gc.T@Hm   # (xu),(yv)
    return KR.reshape(n,n,n,n).transpose(0,2,1,3).reshape(n*n,n*n), KR
def fg(z,Sm,L,P,P2,R,p):
    g,n=L.shape[0],L.shape[1]; s=n*n
    A,B,c,G,H=build(z,L,P,P2,R)
    K,_=Kmat(c,G,H,n); U,S,Vh=np.linalg.svd(Sm-K); f=float(np.sum(S[p:]**2))
    E=((U[:,p:]*(2*S[p:]))@Vh[p:,:]).reshape(n,n,n,n).transpose(0,2,1,3).reshape(s,s)  # (xu),(yv)
    Gm=G.reshape(-1,s); Hm=H.reshape(-1,s); cr=c.reshape(-1)
    gGc=-(E@Hm.T).T            # d f / d (c G)  shape (gR, s)
    gH=-(Gm*cr[:,None])@E      # (gR,s)
    gc=np.sum(gGc*Gm,axis=1).reshape(g,R)
    gG=(gGc*cr[:,None]).reshape(g,R,n,n); gH=gH.reshape(g,R,n,n)
    gA=np.einsum('gik,grij,glj->rkl',L,gG,P,optimize=True); gB=np.einsum('gik,grij,glj->rkl',L,gH,P2,optimize=True)
    return f,np.concatenate([gA.ravel(),gB.ravel(),gc.ravel()])
if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--grp'); ap.add_argument('--mode',default='conj'); ap.add_argument('--mode2',default=None)
    ap.add_argument('--R',type=int,default=1); ap.add_argument('--p',type=int,nargs='+'); ap.add_argument('--seeds',type=int,default=3)
    ap.add_argument('--its',type=int,default=3000); ap.add_argument('--k',type=int,default=1); ap.add_argument('--save',default=None)
    a=ap.parse_args(); T=group(a.grp)
    L,P=perms(T,a.mode,a.k); _,P2=perms(T,a.mode2 or a.mode,a.k); g,n=L.shape[0],L.shape[1]; Sm=swap(n).reshape(n*n,n*n); R=a.R
    for p in a.p:
        for sd in range(a.seeds):
            rng=np.random.default_rng(sd); t=time.time()
            z=np.concatenate([rng.standard_normal(2*R*n*n)/np.sqrt(n),np.ones(g*R)])
            res=minimize(fg,z,args=(Sm,L,P,P2,R,p),jac=True,method='L-BFGS-B',options=dict(maxiter=a.its,maxcor=50,gtol=1e-14,ftol=1e-18))
            A,B,c,G,H=build(res.x,L,P,P2,R); K,KR=Kmat(c,G,H,n)
            sv=np.linalg.svd(Sm-K,compute_uv=False); osr=np.linalg.svd(KR,compute_uv=False); q=int(np.sum(osr>1e-8*osr[0]))
            print(f"{a.grp} k={a.k} N={n} mode={a.mode} R={R} p={p} seed={sd}: tail={res.fun:.3e} sv[p-1],sv[p]={sv[p-1]:.3e},{sv[p]:.3e} gap={sv[p-1]/max(sv[p],1e-300):.2e} q={q} nit={res.nit} t={time.time()-t:.0f}s",flush=True)
            if a.save and res.fun<1e-10: np.save(f"{a.save}_p{p}_s{sd}.npy",res.x)
