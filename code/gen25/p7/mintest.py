import sys, json, time, numpy as np
from spread import *; from load import *
cases = []
for f in ['exact_4_5_3.json', 'exact_5_7_4.json', 'exact_5_6_6.json', 'exact_6_7_7.json']:
    R, C, p, q = load_g19(f); cases.append(('g19 ' + f[6:-5], R, C, p, q))
R, C, p, q = load_g24g(); cases.append(('g24g 6_12_3', R, C, p, q))
R, C, p, q = load_mono('shape_8_8_8_S8T8U0.json'); cases.append(('Q8 8_8_8', R, C, p, q))
R, C = cd(3); cases.append(('CD N=8', R, C, 8, 8))
R, C = cd(4); cases.append(('CD N=16', R, C, 100, 16))
# App B balanced thresholds (h=3 N=3, h=4 N=5)
for h, s0bal in [(3, -1), (4, -3), (5, -3)]:
    N, I, A = appb(h); R = sum(A[s] for s in I if s < s0bal); C = swapT(N) - R
    d = analyze(R, C); cases.append((f'AppB h={h} s0={s0bal}', R, C, d['p'], d['q']))
out = []
for name, R, C, p, q in cases:
    d0 = analyze(R, C, p=p, q=q); t = time.time()
    res = {}
    for nz in [0, 1, 2]:
        best = np.inf
        for seed in range(3 if nz else 2):
            L, R1, C1, h = minimise(R, C, p, q, iters=250, nz=nz, seed=seed, init=(0 if seed == 0 else 0.3))
            # recheck exactly (ranks may grow by nz)
            dd = analyze(R1, C1)
            best = min(best, max(dd['opR'] / np.sqrt(q + nz), dd['opC'] / np.sqrt(p + nz)))
        res[nz] = best
    print(f"{name:18s} p={p} q={q} L0={d0['L']:.4g} Linv={d0['Linv']:.4g} | min over G: {res[0]:.4g}  G+1 move: {res[1]:.4g}  G+2 moves: {res[2]:.4g}  ({time.time()-t:.0f}s)", flush=True)
    out.append(dict(name=name, p=p, q=q, L0=d0['L'], Linv=d0['Linv'], minG=res[0], minG1=res[1], minG2=res[2]))
json.dump(out, open('mintest.json', 'w'), indent=0)
