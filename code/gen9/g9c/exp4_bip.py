"""G9c exp4: the refined benchmark = random bipartite d-regular graph (sides L, R of size N = n/2).
(1) F_t = sum_{T subset R, |T|=t} 2^{|L \\ N(T)|}; i(G) = sum_t F_t; class entropies; which t carry weight.
(2) Costs per profile tau = t d / N (per vertex of one side, leading order, N -> infinity):
      enumeration over T:        (tau/d) log2(e d/tau)
      broadcast network (U/N):   tau log2(e d/tau)        [marks = edges at T, density tau/d]
      Theorem B (class (t, y)):  needs h >= H(e^-tau) N  (dense side) -> Prop. R level
    hybrid: enumerate T1 = T cap R1 (fraction rho of R), network on the rest:
      (tau/d) log2(ed/tau) * [ rho + e_net * d * (1-rho) * exp(-tau*rho) ]
(3) Gottlieb-Kantor rank check: rank of the disjointness matrix of s-subsets vs s'-subsets of [M] is C(M, min(s,s'))
    when s + s' <= M  (so a direct network pays log C(M, min) per cluster pair: no 'O(1) bit per vertex' by change of basis).
"""
import itertools, math, random, sys
import numpy as np, networkx as nx
from math import comb

def bip_regular(N, d, seed):
    rng = random.Random(seed)
    while True:
        perm_ok = True; nb = [set() for _ in range(N)]
        for _ in range(d):
            p = list(range(N)); rng.shuffle(p)
            for u in range(N):
                if p[u] in nb[u]: perm_ok = False
                nb[u].add(p[u])
        if perm_ok: return [sum(1 << y for y in s) for s in nb]   # nb[u] = mask of L-neighbours of u in R

print("(1) exact profile sums, N = 20 per side")
for d in [3, 4, 6]:
    N = 20
    nbm = bip_regular(N, d, 11 + d)
    full = (1 << N) - 1
    cover = np.zeros(1 << N, dtype=np.int64)   # N(T) mask for every T via DP
    for T in range(1, 1 << N):
        b = T & -T
        cover[T] = cover[T ^ b] | nbm[b.bit_length() - 1]
    def popc(a):
        a = a.copy(); c = np.zeros_like(a)
        while a.any():
            c += a & 1; a >>= 1
        return c
    pc_cov = popc(cover); pc_T = popc(np.arange(1 << N, dtype=np.int64))
    Ft = [0] * (N + 1)
    for t in range(N + 1):
        sel = pc_cov[pc_T == t]
        vals, cnts = np.unique(sel, return_counts=True)
        Ft[t] = sum(int(c) * (1 << (N - int(v))) for v, c in zip(vals, cnts))
    tot = sum(Ft)
    print(f" d={d}: log2 i(G) = {math.log2(tot):.2f}  (n = {2*N});  t : log2 C(N,t), log2 F_t, mean |L\\N(T)|")
    for t in range(0, N + 1, 2):
        sel = pc_cov[pc_T == t]
        print(f"    t={t:2d} tau={t*d/N:4.1f}  logC={math.log2(comb(N,t)):6.2f}  logF_t={math.log2(Ft[t]):6.2f}  mean uncovered={N-sel.mean():5.2f}  N e^-tau={N*math.exp(-t*d/N):5.2f}")

print("(2) leading-order costs per side-vertex at tau = 1, and hybrid optimum (units: (1/d) log2(e d))")
for ed in [0.25, 0.5, 1.0, 2.0, 4.0]:
    best = min(((rho + ed * (1 - rho) * math.exp(-rho)), rho) for rho in [i / 1000 for i in range(1001)])
    print(f"   e_net*d = {ed:4.2f}: pure enumeration 1.000, pure network {ed:5.3f}, best hybrid {best[0]:.3f} at rho = {best[1]:.3f}")

print("(3) rank of disjointness of s-sets vs s'-sets of [M]")
for (M, s, sp) in [(8, 2, 3), (8, 1, 4), (9, 2, 4), (9, 3, 3), (10, 2, 5)]:
    A = list(itertools.combinations(range(M), s)); B = list(itertools.combinations(range(M), sp))
    Mx = np.array([[0 if set(a) & set(b) else 1 for b in B] for a in A], dtype=float)
    print(f"   M={M} s={s} s'={sp}: rank = {np.linalg.matrix_rank(Mx)}, C(M,min) = {comb(M, min(s, sp))}")
