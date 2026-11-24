#!/usr/bin/env python3
"""G14e-2: 2 x h grid, lam=1.  T = product of singlets on all diagonal edges (= product of equality tests up
to a local sign flip on one end of each edge, which changes no cut rank).  Candidate family: ALL perfect-matching
(singlet-pairing) states of the 2h(h-1) qubits (they span the SL2-invariant space; product expansions built from the
one-bit identity  delta delta = delta delta - eps eps  are the 2^{#rect} matchings that re-pair inside rectangles).
Costs of a matching M:  r(M) = # pairs crossing the row cut;  c(M) = # pairs joining two different columns."""
import sys, itertools, numpy as np
from scipy.optimize import linprog

h = int(sys.argv[1]) if len(sys.argv) > 1 else 3
nodes = [(a, j) for a in range(2) for j in range(h)]
qubits = [(nd, (1-nd[0], j2)) for nd in nodes for j2 in range(h) if j2 != nd[1]]   # (node, partner)
Q = len(qubits); qi = {q: i for i, q in enumerate(qubits)}
row = [q[0][0] for q in qubits]; col = [q[0][1] for q in qubits]
Z = ((np.arange(2**Q)[:, None] >> np.arange(Q)[None, :]) & 1).astype(np.int8)       # all bit strings
def mstate(pairs):
    v = np.ones(2**Q, dtype=np.int8)
    for p, q in pairs:
        v = v * (Z[:, p] != Z[:, q]) * (1 - 2*Z[:, p])
    return v
def matchings(pts):
    if not pts: yield []; return
    p = pts[0]
    for k in range(1, len(pts)):
        for rest in matchings(pts[1:k] + pts[k+1:]):
            yield [(p, pts[k])] + rest
Tpairs = sorted({tuple(sorted((qi[q], qi[(q[1], q[0])]))) for q in qubits})
T = mstate(Tpairs).astype(float)
ALL = list(matchings(list(range(Q))))
print(f"h={h}: {Q} qubits, {len(ALL)} matchings; T has r={sum(row[p]!=row[q] for p,q in Tpairs)}, c={sum(col[p]!=col[q] for p,q in Tpairs)}")
rc = np.array([(sum(row[p] != row[q] for p, q in m), sum(col[p] != col[q] for p, q in m)) for m in ALL])
rng = np.random.default_rng(1)
# coordinates: project on an orthonormal basis of the span (dimension = Catalan(Q/2))
G = rng.normal(size=(2**Q, 400))
V = np.zeros((len(ALL), 400)); 
for i, m in enumerate(ALL): V[i] = mstate(m).astype(float) @ G
t = T @ G
dim = np.linalg.matrix_rank(V[rng.choice(len(ALL), 2000, replace=False)])
print("dim span(matchings) =", dim)
Uo, s, Wt = np.linalg.svd(V[rng.choice(len(ALL), 3000, replace=False)], full_matrices=False)
Bz = Wt[:dim].T                      # 400 x dim
Vc = V @ Bz; tc = t @ Bz            # exact coordinates in the span
def in_span(mask):
    if mask.sum() == 0: return 1.0
    x, res, rk, sv = np.linalg.lstsq(Vc[mask].T, tc, rcond=None)
    return np.linalg.norm(Vc[mask].T @ x - tc) / np.linalg.norm(tc)
full = rc[:, 0].max()
print("\n(i) max-own-bits feasibility: T in span{M: r(M)<=R} + span{M: c(M)<=C} ?   (rows R=0..,cols C=0..)  1=yes")
for R in range(-1, full+1):
    line = []
    for C in range(-1, full+1):
        mask = (rc[:, 0] <= R) | (rc[:, 1] <= C)
        line.append('1' if in_span(mask) < 1e-8 else '.')
    print(f"   R={R:2d}: " + ' '.join(line))
nrect = h*(h-1)//2
print(f"   product (one-bit identity per rectangle): terms (r,c)=(2i,2({nrect}-i)); feasible iff some i has R>=2i and C>=2({nrect}-i-1)")
# (ii) weighted l1 LP:  min sum_i w_i |x_i|,  w_i = min(2^r, 2^c)
w = np.minimum(2.0**rc[:, 0], 2.0**rc[:, 1])
n = len(ALL)
res = linprog(np.concatenate([w, w]), A_eq=np.hstack([Vc.T, -Vc.T]), b_eq=tc, bounds=(0, None), method='highs')
x = res.x[:n] - res.x[n:]
sup = np.where(np.abs(x) > 1e-7)[0]
print(f"\n(ii) LP (weighted l1): status={res.status}, objective={res.fun:.4f}; support size {len(sup)}; "
      f"true cost sum_supp min(2^r,2^c) = {w[sup].sum():.0f}   [trivial {2.0**full:.0f}]")
for i in sup: print(f"     coef {x[i]:+.4f}  (r,c)={tuple(rc[i])}  pairs={[(qubits[p][0],qubits[q][0]) for p,q in ALL[i]]}")
np.save(f"lp_h{h}_x.npy", x)
# grouped ranks of the LP solution
def ranks(vec):
    Tn = vec.reshape((2,)*Q)   # axis k <-> qubit Q-1-k
    ax = lambda q: Q-1-q
    rq = [ax(q) for q in range(Q) if row[q] == 0]; oq = [ax(q) for q in range(Q) if row[q] == 1]
    rr = np.linalg.matrix_rank(Tn.transpose(rq+oq).reshape(2**len(rq), -1), tol=1e-8)
    cr = []
    for j in range(h):
        cq = [ax(q) for q in range(Q) if col[q] == j]; oq = [ax(q) for q in range(Q) if col[q] != j]
        cr.append(np.linalg.matrix_rank(Tn.transpose(cq+oq).reshape(2**len(cq), -1), tol=1e-8))
    return rr, cr
A = np.zeros(2**Q); B = np.zeros(2**Q)
for i in sup:
    if rc[i, 0] <= rc[i, 1]: A += x[i]*mstate(ALL[i])
    else: B += x[i]*mstate(ALL[i])
print("   check A+B=T:", np.abs(A+B-T).max(), " grouped: row-rank(A)=", ranks(A)[0], " column single-cut ranks of B=", ranks(B)[1],
      " | T: row-rank", ranks(T)[0], "col cut ranks", ranks(T)[1])
