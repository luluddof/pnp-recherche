# min p for P_3 (3-cycle, ring columns) for bond triples up to rotation; numerics (vp.solve)
import sys, itertools, numpy as np, os
from vp import Cyc3, solve
N=int(sys.argv[1]); bmax=int(sys.argv[2]); starts=int(os.environ.get('STARTS','8'))
seen=set()
for b in itertools.product(range(1,bmax+1),repeat=3):
    can=min([b[i:]+b[:i] for i in range(3)])
    if can in seen: continue
    seen.add(can)
    prob=Cyc3(N,can); last=None
    for p in range(N**3,-1,-1):
        f,z=solve(prob,p,starts=starts)
        if f>1e-8: break
        last=(p,f)
    print(f"N={N} bonds={can} q={np.prod(can)} minp={last[0]} (tail {last[1]:.1e}; next p tail {f:.2e})",flush=True)
