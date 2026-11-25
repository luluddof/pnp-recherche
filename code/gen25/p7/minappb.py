import sys, json, time, numpy as np
from spread import *; from load import *
out = []
for h, s0bal in [(2, -1), (3, -1), (4, -3), (5, -3)]:
    N, I, A = appb(h); R = sum(A[s] for s in I if s < s0bal); C = swapT(N) - R
    d0 = analyze(R, C); p, q = d0['p'], d0['q']; t = time.time(); res = {}
    for nz in [0, 1]:
        best = np.inf; bestnorms = None
        for seed in range(2):
            L, R1, C1, hh = minimise(R, C, p, q, iters=400, nz=nz, seed=seed, init=(0 if seed == 0 else 0.3))
            dd = analyze(R1, C1, p=p + nz, q=q + nz)
            if dd['L'] < best: best, bestnorms = dd['L'], (dd['opR'], dd['opC'])
        res[nz] = (best, bestnorms)
    print(f"AppB h={h} N={N} s0={s0bal} p={p} q={q} |R|={d0['opR']:.4g} |C|={d0['opC']:.4g} L0={d0['L']:.4g} Linv={d0['Linv']:.4g} | minG {res[0][0]:.4g} norms {res[0][1][0]:.4g},{res[0][1][1]:.4g} | G+1 {res[1][0]:.4g} ({time.time()-t:.0f}s)", flush=True)
    out.append(dict(h=h, N=N, s0=s0bal, p=p, q=q, L0=d0['L'], Linv=d0['Linv'], minG=res[0][0], minG1=res[1][0]))
json.dump(out, open('minappb.json', 'w'), indent=0)
