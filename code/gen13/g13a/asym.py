from math import comb, log, sqrt, pi, exp
from fractions import Fraction as Fr
def Nm(l,d):
    m=(l-d-1)//2
    if m<0: return 0
    c=0; tot=0
    for j in range(m+1):
        c+=comb(l,j); tot+=c
    return tot
def Nasym(l,d):
    best=0
    for mc in range(0,(l-1-d)//2+1):
        ma=l-1-d-mc
        best=max(best,sum(comb(l,w)*(ma+1-w) for w in range(mc+1)))
    return best
def layered(l,d):
    return min(sum(comb(l-1,c) for c in range(r,l,d+1)) for r in range(d+1))
print("l d | LB=N_m/(l2^l) | asym LB | est sqrt(l)e^{-d^2/2l}*sqrt(2/pi)/(4d^2)*l... | UB layered | sqrt(2/(pi l))e^{-d^2/2l} | UB/LB")
for l in (100,400,1000,4000):
    for g in (1.0,1.5,2.5):
        d=round(sqrt(2*g*l*log(l)))
        lb=Fr(Nm(l,d),l<<l); lb2=Fr(Nasym(l,d),l<<l); ub=Fr(layered(l,d),1<<l)
        m=(l-d-1)//2
        est=comb(l,m)/2**l*(l/(2*d))**2/l
        print(f"{l} {d} | {float(lb):.3e} | {float(lb2):.3e} | C(l,m)2^-l (l/2d)^2/l={est:.3e} | {float(ub):.3e} | {sqrt(2/(pi*l))*exp(-d*d/(2*l)):.3e} | {float(ub/lb):.1f}")
print("threshold d for rho<=l^-C: d_LB (least d with N_m/(l2^l)<=l^-C), d_UB (layered), sqrt((2C-1) l ln l)")
for C2 in (3,4,6):   # 2C
    for l in (100,1000,4000,20000):
        tgt2=Fr(1,l**C2)
        lo,hi=0,l-1
        while lo<hi:
            mid=(lo+hi)//2
            if Fr(Nm(l,mid),l<<l)**2<=tgt2: hi=mid
            else: lo=mid+1
        dl=lo
        du=next(d for d in range(max(1,dl),l) if Fr(layered(l,d),1<<l)**2<=tgt2)
        print(f"C={C2/2} l={l}: d_LB={dl} d_UB={du} ratio={du/dl:.3f} formula={sqrt((C2-1)*l*log(l)):.1f}")
