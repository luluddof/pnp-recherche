# alpha(mu) for P_3 frontiers (B2/U1 formula): alpha(mu)=sup{a: h(mu) <= ||(mu,1-mu)||_{1/(1-a)}}, alpha* = sup_mu.
# Frontier = Pareto boundary; h(mu) = max over boundary points (u,v) of mu u + (1-mu) v.
import numpy as np
from mpmath import mp, mpf, findroot
mp.dps=40
phi=np.linspace(0,np.pi/2,400001)
def frontier(name):
    c,s=np.cos(phi),np.sin(phi)
    if name=='swap_circle': u,v=c,s                           # SWAP App B (normalised u=theta/2,v=T/2)
    elif name=='p3_composite': u,v=(1+2*c)/3,4*s/3            # image of App B under SWAP->P_3 composition
    elif name=='p3_comp_line':                                 # image of the SWAP LINE (if SWAP reached Theorem C)
        t=np.linspace(0,1,400001); u,v=(1+2*t)/3,4*(1-t)/3
    elif name=='p3_bound':                                     # X1 lower bound u >= 1-(2/3)v (best conceivable)
        t=np.linspace(0,1.5,400001); u,v=np.maximum(1-2*t/3,0),t
    P=np.stack([u,v],1)
    P=np.vstack([P,[[0,1],[1,0]]]); P=P[(P[:,1]<=1+1e-12)&(P[:,0]<=1+1e-12)]
    return P
def alpha_mu(h,mu):
    f=lambda a: (mpf(mu)**(1/(1-a))+mpf(1-mu)**(1/(1-a)))**(1-a)
    if f(0)<h: return 0.0
    lo,hi=mpf(0),mpf(1)-mpf(10)**-30
    for _ in range(120):
        m=(lo+hi)/2
        if f(m)>=h: lo=m
        else: hi=m
    return float(lo)
for name in ['swap_circle','p3_composite','p3_comp_line','p3_bound']:
    P=frontier(name); out=[]
    for mu in [0.5,0.6,0.7,0.8,0.9,0.95,0.99,0.999,0.9999,1-1e-5,1-1e-6]:
        h=max(np.max(mu*P[:,0]+(1-mu)*P[:,1]),max(mu,1-mu))
        out.append((mu,round(h,6),round(alpha_mu(mpf(h),mu),4)))
    print(name, out)
