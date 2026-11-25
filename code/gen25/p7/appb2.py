# integer-exact accumulation version (needed at h = 7: float residual)
import sys, json, time, math, numpy as np
import appb_exact as ab
from spread import *
h = int(sys.argv[1]); Nreq = int(sys.argv[2]) if len(sys.argv) > 2 else None
t0 = time.time(); N, mus, I, T = ab.pieces(h, Nreq); d = math.factorial(h); print(f"h={h} N={N} pieces {time.time()-t0:.1f}s", flush=True)
A = {s: np.array(T[s], dtype=object) for s in I}
Fi = np.zeros((N,)*4, dtype=object)
for x in range(N):
    for y in range(N): Fi[x, y, y, x] = d
assert (sum(A.values()) - Fi == 0).all()
print("max |coef|/h!:", max(abs(a).max() for a in A.values()) / d, flush=True)
Ri = np.zeros((N,)*4, dtype=object); res = []
for s0 in I + [I[-1] + 1]:
    Ci = Fi - Ri
    R = np.array(Ri, dtype=float) / d; C = np.array(Ci, dtype=float) / d
    p = ab.rr(Ri.tolist(), N); q = ab.rc(Ci.tolist(), N)
    dd = analyze(R, C, p=p, q=q, check=False); dd['s0'] = s0
    print(f"s0={s0:3d}", fmt(dd), flush=True); res.append({k: float(v) for k, v in dd.items()})
    if s0 in A: Ri = Ri + A[s0]
json.dump(res, open(f'appb_h{h}_N{N}.json', 'w'))
