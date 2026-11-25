import math, numpy as np, b7lib as B, chain_b7 as C
for zname,zf in (('old',lambda lt,x:B.zeta_old(lt,x)),('new',B.zeta_new)):
    best=None
    for x0 in range(3,16):
        K=2**x0; x=float(x0); e=(K-1)/(2*K)
        for ch in [C.chain_greedy(x,e,zf)]+[C.chain_param(x,e,zf,a0,b0) for a0 in np.arange(-4,16.01,1.0) for b0 in (0,.5,1,1.5,2)]:
            mu,s2=C.summarize(ch)
            if best is None or mu>best[0]: best=(mu,s2,x0)
    mu,s2,x0=best
    print(f"brute base (paper level 0, poly space), zeta={zname}: mu={mu:.4f} x0={x0} single c1 a=2 {s2:.3f} a=4 {s2*math.sqrt(2):.3f}; no-saw a=2 {2*math.sqrt(2/mu):.3f}")
