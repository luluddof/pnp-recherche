import sys; sys.path.insert(0,'.')
from appd_explicit import *
print("== legal eps under (D.6), paper schedule, paper network (C=M0/s=%.2f)"%(W_PAPER/(1+ALPHA)))
for lam in (1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12):
    e,dm=eps_legal(lam)
    print(f" lam={lam:.0e}  dmax={dm}  eps_max={e:.6e}  = 1/{1/e:.1f}   (D.5-literal s^(1/1000) needs eps<=1/48000; u>=lam^(4/5) at d=0 needs lam>=eps^-5)")
print("\n== paper as written: rows")
for nm,kw in [("A  paper as written (eps=2^-17, kappa 24, nu 6, trees+subdiv w=45.1, child 3M/32)",{}),
              ("A' same, b4 model numbers (w=11.78,kappa=18)",dict(w=11.78,kappa_mode='loss6')),]:
    for tau,tn in ((1,'2^n'),(0.386,'PPSZ')):
        th=threshold(lambda L: paper(L,**kw)[0]/1.0, tau)
        print(f" {nm:90s} vs {tn:5s}: lambda* = {th:.3g}")
