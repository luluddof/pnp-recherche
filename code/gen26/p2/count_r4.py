# counting version: time/memory exponents (per lam) of outer-square+sweep and of hybrid recursion, large K
import numpy as np
def sq_counts(K, d):
    # nodes 0..K-1, S = {(u,v): m-d <= u < m <= v < m+d}; sweep order 0..K-1. Return (|S|+max step, max mem)
    m = K//2; best_t = 0; best_m = 0
    for j in range(K):   # processing node j: cut before = {u<j<=v}, new = {(j,v): v>j}; mem after = {u<=j<v}
        # count cut edges u<j<=v not in S; S-edges with u<j<=v: u in [m-d, min(j,m)), v in [max(j,m), m+d)
        def inS(ulo, uhi, vlo, vhi):  # number of S-edges with u in [ulo,uhi), v in [vlo,vhi)
            a = max(0, min(uhi, m) - max(ulo, m-d)); b = max(0, min(vhi, m+d) - max(vlo, m)); return a*b
        cut = j*(K-j) - inS(0, j, j, K)
        new = (K-1-j) - inS(j, j+1, j+1, K)
        mem = (j+1)*(K-j-1) - inS(0, j+1, j+1, K)
        best_t = max(best_t, cut + new); best_m = max(best_m, mem)
    return d*d + best_t, best_m
def dcK(k, memo={}):
    if k <= 1: return 0
    if k not in memo: memo[k] = (k//2)*(k-k//2) + dcK(k//2)
    return memo[k]
for K in (64, 256, 1024):
    print(f"K={K}: sweep s(K)/K^2={(((K+1)**2)//4-1)/K**2:.4f}; D&C d(K)/K^2={dcK(K)/K**2:.4f}")
    worst = -1; pts = []
    for d in range(0, K//2+1, max(1, K//32)):
        t, mm = sq_counts(K, d); T, M = t/K**2, mm/K**2
        worst = max(worst, T - (1-M)/3); pts.append((round(M, 4), round(T, 4)))
    # hybrid: middle cut at top (K^2/4), then square scheme in halves (time adds, mem of half)
    hy = []
    for d in range(0, K//4+1, max(1, K//64)):
        t, mm = sq_counts(K//2, d); T = (K*K/4 + t)/K**2; M = mm/K**2
        worst = max(worst, T - (1-M)/3); hy.append((round(M, 4), round(T, 4)))
    print("  square pts (M,T):", pts[::4], "\n  hybrid pts:", hy[::4], f"\n  max (T - (1-M)/3) = {worst:.4f}  (1/(2K)={1/(2*K):.4f})")
