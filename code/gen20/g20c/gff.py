# Profile of free vertex algebras: g bosonic free fields of weight w; character prod_{n>=w}(1-q^n)^{-g}.
# Framework cost (App. B type): cost(h) = max_s [log D(2h+s) + log D(2h-s)] / log dim V_h, D(E)=sum_{e<=E} dim V_e.
# ratio = cost/2 ; App.B (g=w=1) -> sqrt2 asymptotically.  rho_framework = cost/4.
import math, sys
def dims(g,w,Emax):
    sig=[0]*(Emax+1)
    for n in range(w,Emax+1):
        for j in range(n,Emax+1,n): sig[j]+=g*n
    c=[0]*(Emax+1); c[0]=1
    for E in range(1,Emax+1):
        c[E]=sum(sig[j]*c[E-j] for j in range(1,E+1))//E
    return c
def lg(x): return math.log2(x) if x>0 else 0.0
res=[]
for g in [1,2,4,16,64,256,4096]:
  for w in [1,2,3,5,8]:
    Emax=4*120
    c=dims(g,w,Emax)
    cum=[];t=0
    for x in c: t+=x; cum.append(t)
    best=None
    for h in range(max(2,w),121):
        if c[h]<2: continue
        Lh=lg(c[h])
        cost=max(lg(cum[2*h+s])+lg(cum[2*h-s]) for s in range(-2*h,2*h+1))
        r=cost/(2*Lh)
        if Lh>=10 and (best is None or r<best[0]): best=(r,h,Lh)
    print(g,w,best,flush=True)
