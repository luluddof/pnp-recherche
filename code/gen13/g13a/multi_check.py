"""multi-point version over GF(p): points tau_k with directions c(k) (arbitrary multiplicities, some directions unused).
A(ma,mc)={a_x=prod_{k:x_{c(k)}=1}(1-t/tau_k) r_x, deg<=ma, support mu(x)<=mc}, ma+mc+d<=N-1.  Check injectivity of
a -> (a(tau_k) J(tau_k)[x_{c(k)}=0, y_{c(k)}=1])_k and dim = sum_{mu(x)<=mc}(ma+1-mu(x)) (ma=mc=m: sum_x max(0,m+1-mu(x)))."""
import random
from exact_check import rank_mod, polymul, peval, Jeval, randomJ
def run(l,d,cdir,p,ma,mc):
    N=len(cdir); n=1<<l; tau=random.sample(range(1,p),N); Bc=randomJ(l,d,p)
    B=[]
    for x in range(n):
        ks=[k for k in range(N) if (x>>cdir[k])&1]; mu=len(ks)
        if mu>mc or mu>ma: continue
        pi=[1]
        for k in ks: pi=polymul(pi,[1,(-pow(tau[k],p-2,p))%p],p)
        for i in range(ma-mu+1): B.append((x,[0]*i+pi))
    rows=[[] for _ in B]
    for k in range(N):
        c=cdir[k]; Jt=Jeval(Bc,tau[k],p,n); ys=[y for y in range(n) if (y>>c)&1]
        for bi,(x,pol) in enumerate(B):
            if (x>>c)&1: assert peval(pol,tau[k],p)==0; rows[bi]+=[0]*len(ys)
            else:
                v=peval(pol,tau[k],p); rows[bi]+=[v*Jt[x][y]%p for y in ys]
    return len(B),rank_mod(rows,p)
random.seed(5); ok=True; cnt=0
for p in (13,101,(1<<31)-1):
    for l,cdir in ((2,[0,0,1]),(2,[0,0,0,1]),(3,[0,0,1,2]),(3,[0,1,1,1,2]),(3,[0,0,1,1]),(3,[0,1]),(3,[0,0,0,1,1,2,2]),(2,[0,0,1,1,1,1])):
        N=len(cdir)
        for d in range(0,N):
            for mc in range(0,N):
                ma=N-1-d-mc
                if ma<mc: break
                dim,rk=run(l,d,cdir,p,ma,mc); cnt+=1
                if dim!=rk: ok=False; print("VIOLATION",p,l,cdir,d,ma,mc,dim,rk)
print("multi-point injectivity over GF(p), all cases:",ok,cnt,"instances")
