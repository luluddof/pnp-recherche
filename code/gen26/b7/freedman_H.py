# [H] what a Freedman (variance-based) version of R2 would give: zeta_H = 4/t + sqrt(L2)/t + (2/3)(4 sqrt(tL/2)+1) L2/t^2 + 1/(K-1)
import math, numpy as np, b7lib as B, chain_b7 as C
LN2=math.log(2)
def zeta_H(lt, x):
    t=2.0**lt; L=math.log(4)+2*(lt+x)*LN2; L2=2*x*LN2
    return 4/t + math.sqrt(L2)/t + (2/3)*(4*math.sqrt(t*L/2)+1)*L2/t**2 + (1/(2.0**x-1) if x<1000 else 0)
for label, zf in (("R2 (proved)", B.zeta_new), ("R2-Freedman [H]", zeta_H), ("zeta = 4/t only [H, floor of this route]", lambda lt,x: 4/2.0**lt)):
    best=None
    for name,x,e in C.starts('row-R3'):
        for ch in [C.chain_greedy(x,e,zf)]+[C.chain_param(x,e,zf,a0,b0) for a0 in np.arange(-4,16.01,1.0) for b0 in (0,.5,1,1.5,2)]:
            mu,s2=C.summarize(ch)
            if best is None or mu>best[0]: best=(mu,s2,name)
    mu,s2,name=best
    print(f"row-R3 + {label:42s}: mu={mu:.4f} ({name}) single c1 a=2 {s2:.3f}; no-saw a=2 {2*math.sqrt(2/mu):.3f}", flush=True)
