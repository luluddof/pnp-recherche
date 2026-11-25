# Column relaxation of the Bell model: d_c = dist_Ham(chi_c, span_R chi_T) computed EXACTLY by circuit enumeration:
# min-support vectors of span{chi_b : b in T u {c}} with nonzero c-coefficient are determined by |T| independent zeros.
# usage: python3 rcol.py k q ntrials
import numpy as np, itertools, sys
k, q, ntr = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); N = 2**k
rng = np.random.default_rng(7)
pts = np.array(list(itertools.product([0,1], repeat=k)))
H = (-1.0)**(pts @ pts.T % 2)                      # H[x, b] = chi_b(x)
Zsets = np.array(list(itertools.combinations(range(N), q)))   # candidate zero sets (size q)
def dist(T, c):
    B = list(T) + [c]
    M = H[:, B]                                     # N x (q+1)
    sub = M[Zsets]                                  # (#Z, q, q+1)
    # kernel vector of each q x (q+1) block via SVD (last right singular vector)
    _, s, Vt = np.linalg.svd(sub)
    v = Vt[:, -1, :]                                # (#Z, q+1)
    ok = (s[:, -1] > 1e-9)                          # rank q (unique kernel)
    vecs = v @ M.T                                  # (#Z, N)
    w = (np.abs(vecs) > 1e-9).sum(axis=1)
    good = ok & (np.abs(v[:, -1]) > 1e-9)
    return int(w[good].min()) if good.any() else N
def total(T):
    Ts = set(T)
    return sum(dist(T, c) for c in range(N) if c not in Ts)
def idx(vec): return int("".join(map(str, vec)), 2)
res = []
# structured: Hamming ball pieces, subspace, random
cands = []
ball = sorted(range(N), key=lambda b: bin(b).count('1'))
cands.append(("ball-prefix", ball[:q]))
cands.append(("subspace-ish", list(range(q))))
for t in range(ntr):
    cands.append(("random", list(rng.choice(N, size=q, replace=False))))
best = None
for name, T in cands:
    tot = total(T)
    res.append((tot, name, T))
    if best is None or tot < best[0]: best = (tot, name, T)
res.sort()
print(f"k={k} N={N} q={q}: line N^2/(4q) = {N*N/(4*q):.1f}, (N-q)N/(q+1) = {(N-q)*N/(q+1):.1f}")
for tot, name, T in res[:5]: print(f"   sum_c d_c = {tot}  ({name}) T={sorted(int(t) for t in T)}")
