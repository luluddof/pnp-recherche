import sys, numpy as np
from vp import Cyc3, solve
N=int(sys.argv[1])
for spec in sys.argv[2:]:
    b,start=spec.split(':'); b=tuple(map(int,b.split(','))); prob=Cyc3(N,b); last=None
    for p in range(int(start),-1,-1):
        f,z=solve(prob,p,starts=12,seed=p)
        if f>1e-7: break
        last=(p,f)
    print(f"N={N} bonds={b} numerical min p = {last[0] if last else '>'+start} (tail {last[1] if last else 0:.1e}); p-1 tail {f:.2e}",flush=True)
