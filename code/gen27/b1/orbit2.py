# orbit ansatz with per-term diagonal twists: G_g = diag(d_g) L_g A P_g diag(e_g), H_g = diag(d'_g) L_g B P'_g diag(e'_g)  (R seeds)
import numpy as np, argparse, time
from scipy.optimize import minimize
from orbit import group, perms, swap
def unpack(z,g,n,R):
    o=0; A=z[o:o+R*n*n].reshape(R,n,n); o+=R*n*n; B=z[o:o+R*n*n].reshape(R,n,n); o+=R*n*n
    d=z[o:o+g*R*n].reshape(g,R,n); o+=g*R*n; e=z[o:o+g*R*n].reshape(g,R,n); o+=g*R*n
    d2=z[o:o+g*R*n].reshape(g,R,n); o+=g*R*n; e2=z[o:o+g*R*n].reshape(g,R,n)
    return A,B,d,e,d2,e2
def build(z,L,P,P2,R):
    g,n=L.shape[0],L.shape[1]; A,B,d,e,d2,e2=unpack(z,g,n,R)
    M=np.einsum('gik,rkl,glj->grij',L,A,P); M2=np.einsum('gik,rkl,glj->grij',L,B,P2)
    G=d[:,:,:,None]*M*e[:,:,None,:]; H=d2[:,:,:,None]*M2*e2[:,:,None,:]
    return A,B,d,e,d2,e2,M,M2,G,H
def fg(z,F,L,P,P2,R,p):
    g,n=L.shape[0],L.shape[1]; s=n*n
    A,B,d,e,d2,e2,M,M2,G,H=build(z,L,P,P2,R)
    K=np.einsum('grxu,gryv->xyuv',G,H); Mx=(F-K).reshape(s,s)
    U,S,Vh=np.linalg.svd(Mx); f=float(np.sum(S[p:]**2))
    E=((U[:,p:]*(2*S[p:]))@Vh[p:,:]).reshape(n,n,n,n)
    gG=-np.einsum('xyuv,gryv->grxu',E,H); gH=-np.einsum('xyuv,grxu->gryv',E,G)
    gd=np.einsum('grij,grij,grj->gri',gG,M,e); ge=np.einsum('grij,gri,grij->grj',gG,d,M)
    gd2=np.einsum('grij,grij,grj->gri',gH,M2,e2); ge2=np.einsum('grij,gri,grij->grj',gH,d2,M2)
    gM=d[:,:,:,None]*gG*e[:,:,None,:]; gM2=d2[:,:,:,None]*gH*e2[:,:,None,:]
    gA=np.einsum('gik,grij,glj->rkl',L,gM,P); gB=np.einsum('gik,grij,glj->rkl',L,gM2,P2)
    return f,np.concatenate([gA.ravel(),gB.ravel(),gd.ravel(),ge.ravel(),gd2.ravel(),ge2.ravel()])
if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--grp'); ap.add_argument('--mode',default='conj'); ap.add_argument('--mode2',default=None)
    ap.add_argument('--R',type=int,default=1); ap.add_argument('--p',type=int,nargs='+'); ap.add_argument('--seeds',type=int,default=3)
    ap.add_argument('--its',type=int,default=3000); ap.add_argument('--k',type=int,default=1); ap.add_argument('--save',default=None)
    a=ap.parse_args(); T=group(a.grp)
    L,P=perms(T,a.mode,a.k); _,P2=perms(T,a.mode2 or a.mode,a.k); g,n=L.shape[0],L.shape[1]; F=swap(n); R=a.R
    for p in a.p:
        for sd in range(a.seeds):
            rng=np.random.default_rng(sd); t=time.time()
            z=np.concatenate([rng.standard_normal(2*R*n*n)/np.sqrt(n),1+0.3*rng.standard_normal(4*g*R*n)])
            res=minimize(fg,z,args=(F,L,P,P2,R,p),jac=True,method='L-BFGS-B',options=dict(maxiter=a.its,maxcor=50,gtol=1e-14,ftol=1e-18))
            A,B,d,e,d2,e2,M,M2,G,H=build(res.x,L,P,P2,R)
            sv=np.linalg.svd((F-np.einsum('grxu,gryv->xyuv',G,H)).reshape(n*n,n*n),compute_uv=False)
            osr=np.linalg.svd(np.einsum('grxu,gryv->xuyv',G,H).reshape(n*n,n*n),compute_uv=False); q=int(np.sum(osr>1e-8*osr[0]))
            print(f"{a.grp} k={a.k} N={n} mode={a.mode} R={R} p={p} seed={sd}: tail={res.fun:.3e} sv[p-1],sv[p]={sv[p-1]:.3e},{sv[p]:.3e} gap={sv[p-1]/max(sv[p],1e-300):.2e} osr={q} nit={res.nit} t={time.time()-t:.0f}s",flush=True)
            if a.save and res.fun<1e-12: np.save(f"{a.save}_p{p}_s{sd}.npy",res.x)
