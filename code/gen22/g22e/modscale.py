# scaling families (V^{(x)n}; M^{(x)n}), per-copy exact reduction:
#  4rho >= min_{u,v>=0} 2 s_M(u+v) / min( s_M(u), s_V(v), s_M(2u+2v)/2 )      (K = 0 vacuum normalisation)
import math, numpy as np
from chars import *
from scaling import mk_phi
def sl2(N,k,j2):   # bottom-normalised dims of L_k(j), j2 = 2j
    num=[0]*(N+1)
    for n in range(-100,101):
        e=(k+2)*n*n+(j2+1)*n
        if 0<=e<=N: num[e]+=(j2+1+2*n*(k+2))
    return mul(num,heis(N,3),N)
def entropy(d,den=1,mult0=1.0):
    phi,mean=mk_phi(d,den)
    bs=np.exp(np.linspace(math.log(40),math.log(0.02),2500))
    B=[b for b in bs if mean(b)[1]<1e-12]
    E=np.array([mean(b)[0] for b in B]); S=np.array([phi(b)+math.log(mult0)+b*e for b,e in zip(B,E)])
    o=np.argsort(E); E,S=E[o],S[o]
    s0=math.log(d[0]*mult0)
    f=lambda e: np.interp(e,np.concatenate([[0],E]),np.concatenate([[s0],S]))
    return f,E.max()
def best(name,dV,dM,den=1,multV=1.0,multM=1.0):
    sV,eV=entropy(dV,den,multV); sM,eM=entropy(dM,den,multM)
    emax=min(eV,eM)/2.2
    us=np.concatenate([[0],np.exp(np.linspace(math.log(1e-3),math.log(emax/2),120))])
    bestv=(9,0,0)
    for u in us:
        for v in us:
            if u+v>emax: continue
            lam=min(sM(u),sV(v),sM(2*u+2*v)/2)
            if lam<=1e-9: continue
            r=2*sM(u+v)/lam
            if r<bestv[0]: bestv=(r,u,v)
    print("%-36s min 4rho = %.4f at (u,v)=(%.3g,%.3g)  [2sqrt2=2.8284]  s_M(0)=%.3f"%(name,bestv[0],bestv[1],bestv[2],sM(0)))
if __name__=="__main__":
  N=300
  for k in [1,2,4,8]:
      V=sl2(N,k,0)
      for j2 in range(0,k+1):
          best("sl2 level %d, M=L(%g)"%(k,j2/2),V,sl2(N,k,j2))
  nsf=NSferm(2*N); iv=[nsf[2*k] for k in range(N+1)]
  best("Ising, M=sigma",iv,ramond(N)); best("Ising, M=eps",iv,[nsf[2*k+1] for k in range(N)])
  best("Yang-Lee, M=phi(-1/5)",RR(N,{2,3}),RR(N,{1,4}))
  best("NS fermion, M=Ramond(x sqrt2)",nsf,[x for i,x in enumerate(ramond(2*N)) for _ in ([0] if True else [])] if False else None) if False else None
  # NS (half units) vs Ramond: put Ramond in half units
  R=ramond(N); Rh=[0]*(2*N+1)
  for i,x in enumerate(R): Rh[2*i]=x
  best("NS fermion, M=Ramond (x sqrt2)",nsf,Rh,2,1.0,math.sqrt(2))
  best("Heisenberg, M=V",heis(N),heis(N))
