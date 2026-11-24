# Hybrid H_K: Kronecker of App. B pieces A_s (h=4, N1=5) with Hrubes's T1, T2 (n=8), graded by two thresholds:
#   R = sum_{s<s1} A_s (x) T1 + sum_{s<s2} A_s (x) T2,  C = (h! * 2den) X_40 - R.   exact ranks (ev.rank, mod-p for big blocks)
import sys, numpy as np, scipy.sparse as sp, json, time
sys.path.insert(0, '/home/user/pnp-recherche/code/gen20/g20b')
import appb_exact, ev, hr, pieces
h = 4
N1, mus, I, T = appb_exact.pieces(h)
A = {s: np.array(T[s], dtype=np.int64) for s in I}
A = {s: a for s, a in A.items() if a.any()}
S = hr.setup(3, 2); n = S['n']
P = [pieces.piece_tensors(S, c) for c in range(len(S['Cs']))]
T1 = sum(a for a, _ in P); T2 = sum(b for _, b in P)
N = N1 * n
def kron4(a, b):
    return np.einsum('abcd,efgh->aebfcgdh', a, b).reshape(N, N, N, N)
K1 = {s: kron4(a, T1) for s, a in A.items()}; K2 = {s: kron4(a, T2) for s, a in A.items()}
X = sum(K1.values()) + sum(K2.values())
d = X[0, 0, 0, 0]
chk = np.zeros_like(X)
for x in range(N):
    for y in range(N): chk[x, y, y, x] = d
print("pieces sum to", d, "* X_40 exactly:", (chk == X).all(), "; App-B s with nonzero piece:", sorted(A))
rr = lambda R: ev.rank(sp.csr_matrix(R.reshape(N*N, N*N)))[0]
rc = lambda C: ev.rank(sp.csr_matrix(C.transpose(0, 2, 1, 3).reshape(N*N, N*N)))[0]
ss = sorted(A) + [max(A) + 1]
res = []
t0 = time.time()
for s1 in ss:
    R1 = sum((K1[s] for s in A if s < s1), np.zeros_like(X))
    for s2 in ss:
        R = R1 + sum((K2[s] for s in A if s < s2), np.zeros_like(X))
        res.append((s1, s2, rr(R), rc(X - R)))
    print(s1, [(a, b) for _, _, a, b in res[-len(ss):]], f"{time.time()-t0:.0f}s", flush=True)
json.dump(res, open("hybrid_exact_N40.json", "w"))
pts = sorted(set((a, b) for _, _, a, b in res))
front = [pt for pt in pts if not any(o[0] <= pt[0] and o[1] <= pt[1] and o != pt for o in pts)]
print("Pareto frontier (p,q) at N=40:", front)
