import sys, math; sys.path.insert(0,'.')
from run3 import optG
from model2 import W_INC
def thr(tau, **kw):
    lam=1e6
    for it in range(5):
        best=min(optG(lam,D,r=64,**kw) for D in (2,3))
        new=(best*math.sqrt(lam)/tau)**2
        if abs(new/lam-1)<0.003: lam=new; break
        lam=new
    return lam, best*math.sqrt(lam)/W_INC
rows=[("G  + a_d=a_{d-1}+u_d+c_d, u_d ∝ M_d^{2/3}, depth 2-3, r=64 [rigorous]",{}),
      ("H  + Smith/Mahler modes (nu/4) [conditional: Smith profile conj.]",dict(smith=True)),
      ("I  + kappa=3 (black-box recovery, toy numerics only) [open]",dict(smith=True,kappa=3.0)),
      ("I' + kappa=2 (proved witness floor ~2) [open]",dict(smith=True,kappa=2.0))]
for nm,kw in rows:
    a,c=thr(1.0,**kw); b,_=thr(0.386,**kw)
    print(f"{nm:72s} vs 2^n: {a:.3g}  vs PPSZ: {b:.3g}   E*sqrt(lam)/M0 = {c:.2f}", flush=True)
