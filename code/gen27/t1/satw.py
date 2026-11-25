# weight at which each mode-span block S_r (rho: V_h->V_{h+r}) / S'_r (kap: V_{h-r}->V_h) becomes the full parity Hom
import sys; sys.dont_write_bytecode = True
import numpy as np
from sscan import Setup, rowbasis
from appb_exact import partitions
h = int(sys.argv[1]); Wm = int(sys.argv[2]); Dmax = int(sys.argv[3]); P = 16777213
st = Setup(h, Wm, Wm, P); n = st.n
keys = [('rho', r) for r in st.rr] + [('kap', r) for r in st.rk]
S = {k: np.zeros((0, len(st.hmask[k])), dtype=np.int64) for k in keys}; sat = {}; dims = {k: [] for k in keys}
for w in range(Dmax + 1):
    for k in keys:
        if k in sat: continue
        side, r = k
        new = [st.O.mat(mu, r, h).ravel() if side == 'rho' else st.O.mat(mu, r, h - r).ravel() for mu in partitions(w)]
        S[k] = rowbasis(np.concatenate([S[k], np.array(new, dtype=np.int64) % P]), P); dims[k].append(S[k].shape[0])
        if S[k].shape[0] == len(st.hmask[k]): sat[k] = w
for k in keys:
    side, r = k; src, tgt = (h, h + r) if side == 'rho' else (h - r, h)
    print(f"h={h} {side} r={r:>3} Hom(V_{src},V_{tgt}) dim {len(st.hmask[k]):>5}: saturates at wt {sat.get(k, '>'+str(Dmax))}  (src+tgt={src+tgt})  dims {dims[k]}")
