# fast grid version: per-copy rates for A1^n, xi=(alpha/2)^n (MSHIFT=0: labels in V'; 0.5: labels in V_{alpha/2+L})
# plus option CONF2: (A1+A1)^n with xi=(alpha1/2,0): second copy conformal (no flattening loss there).
import math, sys, numpy as np
MSHIFT=float(sys.argv[1]) if len(sys.argv)>1 else 0.0
MODE=sys.argv[2] if len(sys.argv)>2 else "A1"
K=np.arange(1,20000); S=np.arange(-30,32); Q=np.arange(-62,64)
bs=np.exp(np.linspace(math.log(30),math.log(0.004),500))
def lphi(b): return -np.log1p(-np.exp(-b*K)).sum()
LP=np.array([lphi(b) for b in bs])
def w(s,sh=0.0): return (s+sh)**2-(s+sh)
def wconf(s): return s**2          # conformal A1 sector weight
lZD=np.array([math.log(np.exp(-b*w(S)).sum()) for b in bs])+LP
lZM=np.array([math.log(np.exp(-b*w(S,MSHIFT)-(-b*w(S,MSHIFT)).max()*0).sum()) for b in bs])+LP
lZC=np.array([math.log(np.exp(-b*wconf(S)).sum()) for b in bs])+LP
# shift module weights so that min weight is 0
mmin=w(S,MSHIFT).min()
lZM=np.array([math.log(np.exp(-b*(w(S,MSHIFT)-mmin)).sum()) for b in bs])+LP
# logBC on grid (bx index i for labels, by index j for pairs)
WS=w(S); 
LBC=np.zeros((len(bs),len(bs)))
for i,bx in enumerate(bs):
    EM=np.exp(-bx*(w(Q[:,None]+S[None,:],MSHIFT)-mmin))       # [q, s]: label sector q+s
    EMm=np.exp(-bx*(w(Q[:,None]-1+S[None,:],MSHIFT)-mmin))
    for j,by in enumerate(bs):
        ev=np.exp(-by*WS)
        R=EM@ev; Rm=EMm@ev
        LBC[i,j]=math.log(np.sqrt(R*Rm).sum())+LP[i]+LP[j]
def leg1(lZ,x): return float(np.min(lZ+bs*x))
def flat(u,v): return float(np.min(LBC+bs[:,None]*u+bs[None,:]*v))
if MODE=="A1":
    sCost=lambda x: leg1(lZM,x); sLab=lambda x: leg1(lZM,x); sPair=lambda x: leg1(lZD,x)
    fl=lambda u,v: flat(u,v)/2
else:  # (A1+A1) block, copy 2 conformal: profiles are products, flattening loss only from copy 1
    lZMb=lZM[:,None]+lZC[None,:]  # not used directly
    # 2-variable Legendre over (b1,b2) with one total weight: block profile phi(b) = lZM(b)+lZC(b) (same b, total weight)
    lZblkM=lZM+lZC; lZblkD=lZD+lZC
    sCost=lambda x: leg1(lZblkM,x); sLab=sCost; sPair=lambda x: leg1(lZblkD,x)
    # flattening: logBC_1(bx,by) + lZC(bx) + lZC(by)   (conformal copy: rows=cols, rank <= dims)
    LBC2=LBC+lZC[:,None]+lZC[None,:]
    fl=lambda u,v: float(np.min(LBC2+bs[:,None]*u+bs[None,:]*v))/2
grid=np.exp(np.linspace(math.log(0.03),math.log(40.0),70))
bestL=(9,); bestF=(9,)
for u in grid:
    for v in grid:
        cost=2*sCost(u+v); loose=min(sLab(u),sPair(v),sCost(2*u+2*v)/2); f=fl(u,v)
        r1=cost/loose; r2=cost/min(loose,f)
        if r1<bestL[0]: bestL=(r1,u,v)
        if r2<bestF[0]: bestF=(r2,u,v,loose,f)
print("MODE",MODE,"MSHIFT",MSHIFT)
print("  loose            min 4rho = %.4f at u=%.3f v=%.3f"%bestL)
print("  loose+flattening min 4rho = %.4f at u=%.3f v=%.3f  (loose rate %.4f, flat rate %.4f)"%bestF)
# show the flattening rate vs loose along the diagonal near the loose optimum
for u in [0.3,0.6,1.0,1.11,1.5,2.0,3.0]:
    print("   u=v=%.2f: lab %.4f pair %.4f flat %.4f"%(u,sLab(u),sPair(u),fl(u,u)))
