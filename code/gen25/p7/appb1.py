import sys, json, time, numpy as np
from spread import *; from load import *
h = int(sys.argv[1]); Nreq = int(sys.argv[2]) if len(sys.argv) > 2 else None
t0 = time.time(); N, I, A = appb(h, Nreq); print(f"h={h} N={N} pieces {time.time()-t0:.1f}s", flush=True)
F = swapT(N); assert np.abs(sum(A.values()) - F).max() < 1e-9
R = np.zeros_like(F); res = []
print("max |coef| of pieces (normalised by h!):", max(np.abs(a).max() for a in A.values()))
for s0 in I + [I[-1] + 1]:
    C = F - R
    d = analyze(R, C, inv=True, check=False)
    d['s0'] = s0; print(f"s0={s0:3d}", fmt(d), flush=True); res.append({k: float(v) for k, v in d.items()})
    if s0 in A: R = R + A[s0]
json.dump(res, open(f'appb_h{h}_N{N}.json', 'w'))
