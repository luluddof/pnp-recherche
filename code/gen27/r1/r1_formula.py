import numpy as np, sys
sys.path.insert(0,'.'); from r1dp import Hull, fam
from mpmath import mp, mpf; mp.dps=50
def amu(H,mu):
    h=mpf(H.h(mu)); mu=mpf(mu); f=lambda a:(mu**(1/(1-a))+(1-mu)**(1/(1-a)))**(1-a)
    if f(0)<h: return 0.0
    lo,hi=mpf(0),mpf(1)-mpf(10)**-30
    for _ in range(200):
        m=(lo+hi)/2
        if f(m)>=h: lo=m
        else: hi=m
    return float(lo)
for name in sys.argv[1:]:
    H=Hull(fam(name))
    mus=[0.5,0.6,0.7,0.8,0.9,0.95,0.99,1-1e-3,1-1e-4,1-1e-6,1-1e-8]
    print(name,' '.join(f"{m:.8g}:{amu(H,m):.3f}" for m in mus))
