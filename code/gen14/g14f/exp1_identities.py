# G14f exp1: identities and ranks for random bipartite d-regular graphs (small N).
import itertools, random, math, sys
import numpy as np
from fractions import Fraction
from math import comb

def rand_bip(N, d, rng):
    # union of d random perfect matchings without multi-edges
    while True:
        adjR = [set() for _ in range(N)]
        ok = True
        for _ in range(d):
            for _try in range(200):
                p = list(range(N)); rng.shuffle(p)
                if all(p[v] not in adjR[v] for v in range(N)):
                    for v in range(N): adjR[v].add(p[v])
                    break
            else:
                ok = False; break
        if ok: return [sorted(s) for s in adjR]   # adjR[v] = L-neighbours of v in R

def popc(x): return bin(x).count("1")

def run(N, d, seed):
    rng = random.Random(seed)
    adjR = rand_bip(N, d, rng)
    maskR = [sum(1 << u for u in adjR[v]) for v in range(N)]          # N(v) subset of L
    maskL = [sum(1 << v for v in range(N) if u in adjR[v]) for u in range(N)]  # N(u) subset of R
    full = (1 << N) - 1
    # NT[T] = N(T) as mask over L
    NT = [0] * (1 << N)
    for T in range(1, 1 << N):
        low = (T & -T).bit_length() - 1
        NT[T] = NT[T & (T - 1)] | maskR[low]
    NA = [0] * (1 << N)
    for A in range(1, 1 << N):
        low = (A & -A).bit_length() - 1
        NA[A] = NA[A & (A - 1)] | maskL[low]
    # (0) brute force #IS
    Z_bf = sum(1 << (N - popc(NT[T])) for T in range(1 << N))
    # (a) one-side formula over L
    Z_a = sum(1 << (N - popc(NA[A])) for A in range(1 << N))
    # (b) clause form: monotone d-CNF Phi_R on L-variables, V(xi)= violated clauses, sum 2^{|V|}
    Z_b = 0
    for A in range(1 << N):
        V = sum(1 for v in range(N) if maskR[v] & A == 0)
        Z_b += 1 << V
    # (c) Fourier with q=d+1 roots of unity (floating check)
    q = d + 1
    w = np.exp(2j * np.pi / q)
    Z_c = 0.0
    for A in range(1 << N):
        prod = 1.0
        for v in range(N):
            k = popc(maskR[v] & A)
            ind = sum(w ** (j * k) for j in range(q)) / q   # [k==0]
            prod *= (1 + ind)
        Z_c += prod.real
    # (d) Holant / even-subgraph form, by classes |T|=t
    Zt = [0] * (N + 1); Ht = [Fraction(0)] * (N + 1)
    for T in range(1 << N):
        t = popc(T)
        Zt[t] += 1 << (N - popc(NT[T]))
        ev = 1
        for u in range(N):
            k = popc(maskL[u] & T)
            ev *= (1 << (k - 1)) if k >= 1 else 1
        Ht[t] += Fraction(ev, 1 << (d * t))
    okd = all(Fraction(Zt[t]) == (1 << N) * Ht[t] for t in range(N + 1))
    print(f"N={N} d={d} seed={seed}: #IS={Z_bf}  one-side={Z_a==Z_bf} clause={Z_b==Z_bf} "
          f"fourier_err={abs(Z_c-Z_bf):.2e} holant_classes={okd}")
    # (e) cut ranks of F(A)=2^{|R\N(A)|} across L=L1|L2 (first h vertices vs rest)
    for h in sorted(set([1, 2, 3, N // 2])):
        if h >= N: continue
        m1 = (1 << h) - 1
        F = np.zeros((1 << h, 1 << (N - h)))
        for A1 in range(1 << h):
            for A2 in range(1 << (N - h)):
                A = A1 | (A2 << h)
                F[A1, A2] = float(1 << (N - popc(NA[A])))
        rk = np.linalg.matrix_rank(F)
        # C[B1,B2] = #{T : L\N(T) = B1 u B2}
        cnt = {}
        for T in range(1 << N):
            B = full & ~NT[T]
            cnt[B] = cnt.get(B, 0) + 1
        C = np.zeros((1 << h, 1 << (N - h)))
        for B, c in cnt.items():
            C[B & m1, B >> h] = c
        rkC = np.linalg.matrix_rank(C)
        tr1 = len({NT[T] & m1 for T in range(1 << N)})
        tr2 = len({NT[T] >> h for T in range(1 << N)})
        print(f"   cut |L1|={h}: rank F = {rk}, rank C = {rkC}, #traces N(T)&L1 = {tr1}, #traces on L2 = {tr2}, 2^h = {1<<h}")
    # (e') T-only function restricted to |T cap R1| = j : rank across (R1, R2) with small cluster R1
    return Z_bf

def gottlieb(m, s, s2):
    # ranks for s-sets vs s2-sets of [m]: disjointness, 2^{|X cap Y|}, and zero-pattern lower bound
    X = list(itertools.combinations(range(m), s)); Y = list(itertools.combinations(range(m), s2))
    D = np.array([[1.0 if not set(x) & set(y) else 0.0 for y in Y] for x in X])
    P = np.array([[2.0 ** len(set(x) & set(y)) for y in Y] for x in X])
    rng = np.random.default_rng(1)
    Wt = D * rng.uniform(0.5, 2.0, size=D.shape)     # random matrix with the same zero pattern
    return np.linalg.matrix_rank(D), np.linalg.matrix_rank(P), np.linalg.matrix_rank(Wt)

if __name__ == "__main__":
    for (N, d, seed) in [(8, 3, 1), (10, 3, 2), (10, 4, 3), (12, 4, 4)]:
        run(N, d, seed)
    print("Gottlieb / zero-pattern ranks: m, s, s' : rank DISJ, rank 2^{|X cap Y|}, rank(random same support); C(m,min), C(s+s',s)")
    for (m, s, s2) in [(8, 2, 2), (8, 2, 4), (9, 2, 4), (9, 3, 3), (10, 2, 5)]:
        a, b, c = gottlieb(m, s, s2)
        print(f"   m={m} s={s} s'={s2}: {a} {b} {c} ; C(m,min)={comb(m,min(s,s2))}  C(s+s',s)={comb(s+s2,s)}")
