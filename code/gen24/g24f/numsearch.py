# numerical search for the TRUE max N: sum_g (A B_g Bm^T) (x) (U O_g Y) = I (x) I
import sys, time, numpy as np
REG=1e-7
from scipy.optimize import least_squares
from build import setup
def prep(n,E0s,ha,hb):
    C,A,Bb,Mb,Bs,Os=setup(n,E0s,ha,hb)
    gs=[g for g in C.S]
    Bn=np.array([Bs[g].astype(float) for g in gs]); On=np.array([Os[g].astype(float) for g in gs])
    # drop channels with zero B or zero O
    keep=[k for k in range(len(gs)) if np.abs(Bn[k]).sum()>0 and np.abs(On[k]).sum()>0]
    return Bn[keep],On[keep],len(A),len(Bb),len(Mb)
def solve(Bn,On,dA,dB,dM,N,trials=20,seed=0,tol=1e-10,maxfev=3000):
    rng=np.random.default_rng(seed)
    G=len(Bn); T=np.einsum('xv,uy->xvuy',np.eye(N),np.eye(N))
    sizes=[N*dA,N*dB,N*dM,dM*N]
    def unpack(z):
        o=0; out=[]
        for s,shape in zip(sizes,[(N,dA),(N,dB),(N,dM),(dM,N)]):
            out.append(z[o:o+s].reshape(shape)); o+=s
        return out
    def fun(z):
        A,Bm,U,Y=unpack(z)
        beta=np.einsum('xi,gij,vj->gxv',A,Bn,Bm); Om=np.einsum('uk,gkl,ly->guy',U,On,Y)
        return np.concatenate([(np.einsum('gxv,guy->xvuy',beta,Om)-T).ravel(),REG*z])
    def jac(z):
        A,Bm,U,Y=unpack(z)
        BB=np.einsum('gij,vj->giv',Bn,Bm); AB=np.einsum('xi,gij->gxj',A,Bn)
        beta=np.einsum('xi,giv->gxv',A,BB); Om=np.einsum('uk,gkl,ly->guy',U,On,Y)
        OY=np.einsum('gkl,ly->gky',On,Y); UO=np.einsum('uk,gkl->gul',U,On)
        I=np.eye(N)
        JA=np.einsum('ab,giv,guy->xvuyai',I,BB,Om).reshape(N**4,-1) if False else np.einsum('xa,giv,guy->xvuyai',I,BB,Om).reshape(N**4,-1)
        JB=np.einsum('gxj,vb,guy->xvuybj',AB,I,Om).reshape(N**4,-1)
        JU=np.einsum('gxv,uc,gky->xvuyck',beta,I,OY).reshape(N**4,-1)
        JY=np.einsum('gxv,guk,yd->xvuykd',beta,UO,I).reshape(N**4,-1)
        return np.vstack([np.hstack([JA,JB,JU,JY]),REG*np.eye(len(z))])
    best=1e9; bestz=None
    for t in range(trials):
        z0=rng.standard_normal(sum(sizes))
        meth='lm'
        r=least_squares(fun,z0,jac=jac,method=meth,xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=maxfev)
        c=np.linalg.norm(r.fun[:N**4])
        if c<best: best=c; bestz=r.x; bestnorm=np.abs(r.x).max()
        if best<tol: break
    solve.norm=bestnorm
    return best,bestz,unpack
if __name__=="__main__":
    n,E0s,ha,hb=map(int,sys.argv[1:5]); Nmax=int(sys.argv[5]); trials=int(sys.argv[6]) if len(sys.argv)>6 else 20
    Bn,On,dA,dB,dM=prep(n,E0s,ha,hb)
    print("n=%d E0s=%d ha=%d hb=%d dims A,B,M=%d,%d,%d channels=%d"%(n,E0s,ha,hb,dA,dB,dM,len(Bn)),flush=True)
    for N in range(1,Nmax+1):
        t=time.time(); best,z,_=solve(Bn,On,dA,dB,dM,N,trials=trials,seed=N)
        print("  N=%d  min residual %.3e  maxabs %.1e (%.1fs)"%(N,best,solve.norm,time.time()-t),flush=True)
        if best>1e-6 and N>1: pass
