"""LB (Theorem, G11b): rho_max >= rho_avg >= N_m/(l 2^l), N_m = sum_{j<=m} C(l,<=j), m=floor((l-d-1)/2).
UB (G7a layered, blocks of d+1 levels, best offset): rho_k = sum_{cuts c} C(l-1,c)/2^l."""
from math import comb, log, sqrt
from fractions import Fraction as Fr
def tables(l):
    row=[comb(l,i) for i in range(l+1)]; cum=[0]*(l+2)
    for i in range(l+1): cum[i+1]=cum[i]+row[i]           # cum[j+1]=C(l,<=j)
    dc=[0]*(l+2)
    for j in range(l+1): dc[j+1]=dc[j]+cum[j+1]            # dc[m+1]=N_m
    r1=[comb(l-1,i) for i in range(l)]
    return dc,r1
def LB(l,d,dc):
    m=(l-d-1)//2
    return Fr(dc[m+1],l<<l) if m>=0 else Fr(0)
def UB(l,d,r1):
    best=None
    for r in range(d+1):
        s=sum(r1[c] for c in range(r,l,d+1))
        if best is None or s<best: best=s
    return Fr(best,1<<l)
print("l d | LB rho_avg | UB rho (layered) | UB/LB | 4d^2/l")
for l in (10,20,50,100,400,1000):
    dc,r1=tables(l)
    for g in (0.5,1.0,1.5):
        d=max(1,round(sqrt(2*g*l*log(l))))
        if d>=l: continue
        a=LB(l,d,dc); b=UB(l,d,r1)
        print(f"{l} {d} (gamma={g}) | {float(a):.3e} | {float(b):.3e} | {float(b/a):.2f} | {4*d*d/l:.2f}")
print("\nthreshold degree for rho <= l^-C:  d_LB = least d not excluded by the theorem, d_UB = least d achieved by layered; sqrt((2C-1) l ln l)")
for C in (1.5,2,3):
    for l in (20,50,100,400,1000,4000):
        dc,r1=tables(l); tgt=Fr(1,1)/Fr(l)**int(2*C)  # compare squares to stay exact for half-integers
        dl=next(d for d in range(0,l) if LB(l,d,dc)**2<=tgt)
        du=next(d for d in range(1,l) if UB(l,d,r1)**2<=tgt)
        print(f"C={C} l={l}: d_LB={dl} d_UB={du} ratio={du/dl:.3f}  sqrt((2C-1) l ln l)={sqrt((2*C-1)*l*log(l)):.1f}")
