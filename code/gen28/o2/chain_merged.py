# Chains of independent exact points at the same N: pieces A_1 = R_1, A_s = R_s - R_{s-1}, A_{m+1} = S - R_m.
# k-th tensor power, type classes = compositions kappa of k into m+1 parts. Exact merged row / col ranks (mod 2 primes) vs
#  'ideal' multinom(kappa) prod d_s^kappa_s (d_s = p_s - p_{s-1}; cols e_s = q_{s-1} - q_s)  and
#  'subadditive' multinom(kappa) prod rank(A_s)^kappa_s  (true piece ranks).  Also spectra / commutation of M_s = A_s S^-1.
import sys, itertools, math, numpy as np
from pieces import *
from merged import rk, kron_list
def comps(k, m):
    for c in itertools.product(range(k+1), repeat=m):
        if sum(c) == k: yield c
def classes_sum(mats, k, kappa):
    A = None
    for w in set(itertools.permutations([s for s, c in enumerate(kappa) for _ in range(c)])):
        T = kron_list([mats[s] for s in w]); A = T if A is None else A + T
    return A
files = sys.argv[2:]; k = int(sys.argv[1])
dec = [load_json(f) for f in files]; N = dec[0][0]; S = swap(N); N2 = N*N
Rs = sorted([d[1][0] for d in dec], key=lambda R: rk(rowflat(R))[0])
pieces = [Rs[0]] + [b - a for a, b in zip(Rs, Rs[1:])] + [S - Rs[-1]]
P = [0] + [rk(rowflat(R))[0] for R in Rs] + [N2]; Qs = [N2] + [rk(colflat(S - R))[0] for R in Rs] + [0]
d = [b - a for a, b in zip(P, P[1:])]; e = [a - b for a, b in zip(Qs, Qs[1:])]
rr = [rk(rowflat(A))[0] for A in pieces]; cc = [rk(colflat(A))[0] for A in pieces]
print(f"N={N} chain (p,q) = {list(zip(P,Qs))}; piece row ranks {rr} (ideal d {d}); piece col ranks {cc} (ideal e {e})")
Sr = rowflat(S).astype(float); Sc = colflat(S).astype(float)
Mr = [rowflat(A) @ np.linalg.inv(Sr) for A in pieces]; Mc = [colflat(A) @ np.linalg.inv(Sc) for A in pieces]
comm = max(np.abs(a @ b - b @ a).max() for a in Mr for b in Mr); commc = max(np.abs(a @ b - b @ a).max() for a in Mc for b in Mc)
print(f"  max |[M_s,M_t]| row {comm:.3g} col {commc:.3g}")
rowm = [rowflat(A) for A in pieces]; colm = [colflat(A) for A in pieces]
for kappa in comps(k, len(pieces)):
    mult = math.factorial(k)//math.prod(math.factorial(c) for c in kappa)
    r = rk(classes_sum(rowm, k, kappa)); c = rk(classes_sum(colm, k, kappa))
    idr = mult*math.prod(x**y for x, y in zip(d, kappa)); idc = mult*math.prod(x**y for x, y in zip(e, kappa))
    sr = mult*math.prod(x**y for x, y in zip(rr, kappa)); sc = mult*math.prod(x**y for x, y in zip(cc, kappa))
    print(f"  kappa={kappa}: row {r} (ideal {idr}, subadd {sr}, mult-free {sr//mult})   col {c} (ideal {idc}, subadd {sc}, mult-free {sc//mult})", flush=True)
