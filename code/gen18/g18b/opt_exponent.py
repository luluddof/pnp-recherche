#!/usr/bin/env python3
"""G18b opt_exponent: numerical optimisation of the cost exponent (natural log, divided by n) of the interpolation
route for several charges of the error part.  Model (from P (2.11), (4.16), hD<=2M<=6n, 12 D^2<=n, D=r^{2Q}):
  E/n = 2 ln2 / r  +  3 (R/Q) ln(2eQ/(3R))  +  ERR,
  V0  P as written      : ERR = 6 r tail ln(2 r^{2Q}+1),            tail = 2 exp(-R^2/(8Q))   (P's own bound)
  V0h same, Hoeffding   : tail = 2 exp(-R^2/(2Q))
  V1  reachable, common : ERR = 12 tail (Q ln 2)                    (rank<=2hc, row norm<=2^{Q-1})
  V1p reachable, per-pt : ERR = 12 rho_pw (Q ln 2)
  V2  mod p~R, common   : ERR = 12 tail ln(2R)
  V3  mod 3, per-point  : ERR = 12 ln3 * rho_pw, rho_pw = mass of block-boundary levels (blocks of length R everywhere; ~1/R for R<<sqrt Q, ~rho_c for R>>sqrt Q)
  V3p mod p~R per-point : ERR = 12 ln(2R) * rho_pw
Continuous relaxation (no integrality); asymptotic comparison only."""
import math
def E(v,L,r,Q,R):
    if R>=Q/2 or R<1: return float('inf')
    t8=2*math.exp(-R*R/(8*Q)); t2=2*math.exp(-R*R/(2*Q)); rc=min(1.0,2*math.sqrt(2/(math.pi*Q))*math.exp(-R*R/(2*Q)))
    # per-point piecewise blocks of length R everywhere: rho = mass of the boundary levels (+-R/2, +-3R/2, ...),
    # Gaussian approximation of Bin(Q-1, ~1/2), sigma = sqrt(Q)/2
    sg=math.sqrt(Q)/2; rp=0.0; j=0
    if R<sg: rp=1.0/R      # Poisson summation: relative error exp(-2 pi^2 sg^2/R^2) < 3e-9
    while R>=sg:
        term=2*math.exp(-((j+0.5)*R)**2/(2*sg*sg))/(sg*math.sqrt(2*math.pi)); rp+=term; j+=1
        if term<1e-18*max(rp,1e-300) or j>200000: break
    rp=min(1.0,rp)
    base=2*math.log(2)/r+3*(R/Q)*math.log(2*math.e*Q/(3*R))
    err={'V0':6*r*t8*(2*Q*math.log(r)+math.log(2)),'V0h':6*r*t2*(2*Q*math.log(r)+math.log(2)),
         'V1':12*t2*Q*math.log(2),'V1p':12*rp*Q*math.log(2),'V2':12*t2*math.log(2*R),
         'V3':12*math.log(3)*rp,'V3p':12*math.log(2*R)*rp}[v]
    return base+err
def opt(v,L):
    best=(float('inf'),)
    lr=0.7
    while lr<min(L/8,200):
        r=math.exp(lr); Q=(L-math.log(12))/(4*lr)
        if Q>4 and math.log(2*Q)+2*lr<=L:
            x=0.0
            while True:
                R=math.exp(x)
                if R>=Q/2: break
                e=E(v,L,r,Q,R)
                if e<best[0]: best=(e,r,Q,R)
                x+=0.03
        lr*=1.02
    return best
print(f"{'log2 n':>8} {'variant':>4} {'E/n':>10} {'r':>10} {'Q':>10} {'R':>9} {'R/sqrtQ':>8} | E/n / (L2^2/sqrtL)  E/n / (L2^1.5/sqrtL)  E/n / (L2/sqrtL)")
for lg in (1e4,1e6,1e9,1e12,1e15,1e20,1e30):
    L=lg*math.log(2); L2=math.log(L)
    for v in ('V0','V0h','V1','V1p','V2','V3p','V3'):
        e,r,Q,R=opt(v,L)
        print(f"{lg:8.0e} {v:>4} {e:10.3e} {r:10.3g} {Q:10.3g} {R:9.3g} {R/math.sqrt(Q):8.2f} | {e/(L2*L2/math.sqrt(L)):10.3f} {e/(L2**1.5/math.sqrt(L)):14.3f} {e/(L2/math.sqrt(L)):14.3f}")
    print(f"{'':8} lower bounds (functional of G14d Thm 4): base n: {L2*L2/math.sqrt(2*L):.3e};  base 3 (this note, leading term ln3/sqrt(2 pi) sqrt(L2/L)): {math.log(3)/math.sqrt(2*math.pi)*math.sqrt(L2/L):.3e}")
