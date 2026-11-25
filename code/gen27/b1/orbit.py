# Orbit (twirl) ansatz: K = sum_{g in Gamma} sum_{r<R} c_{g,r} (L_g A_r P_g) (x) (L_g B_r P'_g),
# Gamma acting regularly on [N] (N = |Gamma|), seeds A_r, B_r dense real (NON-monomial), q <= R*|Gamma|.
# mode conj: P_g = L_g^{-1};  left: P_g = I;  lr: P_g = R_g (right mult).  Objective: tail sum_{i>p} sigma_i(S-K)^2 (VP as gen26 vp11).
import numpy as np, sys, argparse, time
from scipy.optimize import minimize

def closure(gens):
    els=[np.eye(gens[0].shape[0],dtype=complex)]
    key=lambda M: tuple(np.round(M,6).ravel().view(float))
    seen={key(els[0]):0}; i=0
    while i<len(els):
        for g in gens:
            M=els[i]@g; k=key(M)
            if k not in seen: seen[k]=len(els); els.append(M)
        i+=1
    n=len(els); T=np.zeros((n,n),dtype=int)
    for a in range(n):
        for b in range(n): T[a,b]=seen[key(els[a]@els[b])]
    return T
def group(name):
    w=lambda k: np.exp(2j*np.pi/k)
    if name.startswith('Z') and 'x' not in name and name[1:].isdigit():
        k=int(name[1:]); return closure([np.array([[w(k)]])])
    if 'x' in name:
        a,b=name.split('x',1); Ta,Tb=group(a),group(b); na,nb=len(Ta),len(Tb)
        T=np.zeros((na*nb,)*2,dtype=int)
        for i in range(na*nb):
            for j in range(na*nb): T[i,j]=Ta[i//nb,j//nb]*nb+Tb[i%nb,j%nb]
        return T
    if name=='Q8': return closure([np.diag([1j,-1j]),np.array([[0,1],[-1,0]],dtype=complex)])
    if name.startswith('Q') : # generalized quaternion of order 4m
        k=int(name[1:])//2; return closure([np.diag([w(k),1/w(k)]),np.array([[0,1],[-1,0]],dtype=complex)])
    if name.startswith('D'): # dihedral of order 2m  (D8 = order 8)
        m=int(name[1:])//2; c,s=np.cos(2*np.pi/m),np.sin(2*np.pi/m)
        return closure([np.array([[c,-s],[s,c]],dtype=complex),np.array([[1,0],[0,-1]],dtype=complex)])
    if name=='H27': # Heisenberg mod 3
        X=np.roll(np.eye(3),1,axis=0).astype(complex); Z=np.diag([1,w(3),w(3)**2]); return closure([X,Z])
    if name=='SL23': return closure([np.array([[1,1],[0,1]],dtype=complex)%3, ]) # placeholder (not used)
    raise ValueError(name)

def perms(T,mode,k=1):
    L,P=perms1(T,mode)
    if k>1: L=np.stack([np.kron(l,np.eye(k)) for l in L]); P=np.stack([np.kron(l,np.eye(k)) for l in P])
    return L,P
def perms1(T,mode):
    n=len(T); L=np.zeros((n,n,n)); P=np.zeros((n,n,n))
    inv=[int(np.where(T[g]==0)[0][0]) for g in range(n)]
    for g in range(n):
        for k in range(n): L[g,T[g,k],k]=1           # L_g e_k = e_{gk}
        if mode=='conj': P[g]=L[g].T
        elif mode=='left': P[g]=np.eye(n)
        elif mode=='lr':
            for k in range(n): P[g,T[k,g],k]=1       # right mult by g: e_k -> e_{kg} ; as matrix acting on right of A
            P[g]=P[g].T
    return L,P

def build(z,L,P,P2,R):
    n=L.shape[1]; g=L.shape[0]
    A=z[:R*n*n].reshape(R,n,n); B=z[R*n*n:2*R*n*n].reshape(R,n,n); c=z[2*R*n*n:].reshape(g,R)
    G=np.einsum('gik,rkl,glj->grij',L,A,P); H=np.einsum('gik,rkl,glj->grij',L,B,P2)
    return A,B,c,G,H
def fg(z,F,L,P,P2,R,p):
    n=L.shape[1]; s=n*n
    A,B,c,G,H=build(z,L,P,P2,R)
    K=np.einsum('gr,grxu,gryv->xyuv',c,G,H); M=(F-K).reshape(s,s)
    U,S,Vh=np.linalg.svd(M); f=float(np.sum(S[p:]**2))
    E=((U[:,p:]*(2*S[p:]))@Vh[p:,:]).reshape(n,n,n,n)
    gG=-np.einsum('xyuv,gryv->grxu',E,H)   # d f/d (c G)
    gH=-np.einsum('xyuv,grxu->gryv',E,G)
    gc=np.einsum('grxu,grxu->gr',gG,G)
    gGc=gG*c[:,:,None,None]; gHc=gH*c[:,:,None,None]
    gA=np.einsum('gik,grij,glj->rkl',L,gGc,P); gB=np.einsum('gik,grij,glj->rkl',L,gHc,P2)
    return f,np.concatenate([gA.ravel(),gB.ravel(),gc.ravel()])
def swap(N):
    X=np.zeros((N,)*4)
    for x in range(N):
        for y in range(N): X[x,y,y,x]=1
    return X
def run(T,mode,mode2,R,p,seed,its,cfix=False,k=1):
    L,P=perms(T,mode,k); _,P2=perms(T,mode2,k); n=L.shape[1]; F=swap(n)
    rng=np.random.default_rng(seed)
    z=np.concatenate([rng.standard_normal(2*R*n*n)/np.sqrt(n),np.ones(L.shape[0]*R)])
    res=minimize(fg,z,args=(F,L,P,P2,R,p),jac=True,method='L-BFGS-B',options=dict(maxiter=its,maxcor=30,gtol=1e-14,ftol=1e-16))
    z=res.x; A,B,c,G,H=build(z,L,P,P2,R)
    K=np.einsum('gr,grxu,gryv->xyuv',c,G,H); s=n*n
    sv=np.linalg.svd((F-K).reshape(s,s),compute_uv=False)
    osr=np.linalg.svd(np.einsum('gr,grxu,gryv->xuyv',c,G,H).reshape(s,s),compute_uv=False)
    q=int(np.sum(osr>1e-8*osr[0]))
    return res.fun,sv,q,z
if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--grp'); ap.add_argument('--mode',default='conj'); ap.add_argument('--mode2',default=None)
    ap.add_argument('--R',type=int,default=1); ap.add_argument('--p',type=int,nargs='+'); ap.add_argument('--seeds',type=int,default=3)
    ap.add_argument('--its',type=int,default=3000); ap.add_argument('--save',default=None); ap.add_argument('--k',type=int,default=1)
    a=ap.parse_args(); T=group(a.grp); n=len(T)*a.k; m2=a.mode2 or a.mode
    for p in a.p:
        best=None
        for sd in range(a.seeds):
            t=time.time(); f,sv,q,z=run(T,a.mode,m2,a.R,p,sd,a.its,k=a.k)
            gap=sv[p-1]/max(sv[p],1e-300) if p>0 else 0
            print(f"{a.grp} N={n} mode={a.mode}/{m2} R={a.R} p={p} seed={sd}: tail={f:.3e} sv[p-1..p+1]={sv[max(p-1,0)]:.3e},{sv[p]:.3e},{sv[p+1]:.3e} gap={gap:.2e} osr={q} t={time.time()-t:.0f}s",flush=True)
            if best is None or f<best[0]: best=(f,z)
        if a.save: np.save(f"{a.save}_p{p}.npy",best[1])
