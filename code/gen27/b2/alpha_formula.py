# alpha* = sup{alpha: exists mu, h(mu) <= ||(mu,1-mu)||_{1/(1-alpha)}}  (self-similar solutions of the Lemma-U level map)
import numpy as np, sys
sys.argv=[sys.argv[0]]; src=open("dp_alpha2.py").read().split("mus=")[0].replace("np.linspace(0,1,2001)","np.linspace(0,1,200001)"); exec(src)
from mpmath import mp, mpf
def alpha_star(P):
    best=(0,None)
    mus=np.linspace(0.5,0.95,901)
    for mu in mus:
        h=np.max(mu*P[:,0]+(1-mu)*P[:,1])
        lo,hi=0.0,1-1e-9
        f=lambda a: (mu**(1/(1-a))+(1-mu)**(1/(1-a)))**(1-a)   # ||.||_q, q=1/(1-a); decreasing in a
        if f(lo)<h: continue
        for _ in range(60):
            m=(lo+hi)/2
            if f(m)>=h: lo=m
            else: hi=m
        if lo>best[0]: best=(lo,mu)
    return best
for name in ['circle','line','ellp1.5','cornerC4','cornerC10','powC0.75','powC0.6']:
    a,mu=alpha_star(family(name)); print(f"{name:10s} alpha*={a:.4f} at mu={mu}")
