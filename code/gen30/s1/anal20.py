import numpy as np, pickle
from core import *
from sl import trunc
k = 3; D = 64; V = np.diag(beta(k)); f = beta(k)
keep = pickle.load(open('keep_r6_t12.pkl', 'rb')); r = 6
rng = np.random.default_rng(0); L = trunc(rng.standard_normal((D, D)), r)
for it in range(20000):
    C = pinv_trans(ptrans(V - L) * keep); L = trunc(V - C, r)
    if np.linalg.norm(V - C - L) < 1e-13: break
E = V - C; s = np.linalg.svd(E, compute_uv=False); print('sing vals', np.round(s[:9], 6))
c = ptrans(C); sup = np.argwhere(np.abs(c) > 1e-9); print('support size', len(sup))
for a, b in sup:
    if a: print('off-diag Pauli a=', a, 'b=', b, 'coef', c[a, b])
e = np.diag(E); print('diag error support', np.flatnonzero(np.abs(e) > 1e-9), 'values', np.round(e[np.abs(e) > 1e-9], 4))
off = E - np.diag(e); nz = np.argwhere(np.abs(off) > 1e-9); print('#offdiag nonzeros', len(nz))
