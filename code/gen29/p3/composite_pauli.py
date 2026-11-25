# Composite C' = K12' K23' in the M-picture with Majorana Bell windows; exact rank(M - C') mod 2 primes.
import numpy as np, sys, math, json
from flint import nmod_mat
from core import *
from bellwin import majorana_windows
def had(k):
    N = 2**k; H = np.zeros((N, N), dtype=np.int64)
    for i in range(N):
        for j in range(N): H[i, j] = (-1)**pc(i & j)
    return H
def build(k, phiL, which, p):
    # integer matrix N*L*K' mod p
    N = 2**k; n = N*N; msk = N - 1; H = had(k); A = np.zeros((n, n), dtype=np.int64)
    ph = np.array([int(v) % p for v in phiL], dtype=np.int64).reshape(N, N)  # ph[hi, lo]
    for y in range(N):
        for x in range(N):
            if which == 12:   # block for fixed y: x' <- x : sum_xi H[xi,x^x'] phi(xi | y<<k)
                col = (H[:, :] * 0)
            pass
    if which == 12:
        for y in range(N):
            D = ph[y, :]               # phi(xi, y) over xi (lo)
            B = (H * D[None, :]) @ H % p   # B[x', x] = sum_xi H[x',xi] D[xi] H[xi,x]
            for xp in range(N):
                A[xp | (y << k), np.arange(N) | (y << k)] = B[xp, :]
    else:
        for x in range(N):
            D = ph[:, x]               # phi(x | eta<<k) over eta (hi)
            B = (H * D[None, :]) @ H % p
            for yp in range(N):
                A[x | (yp << k), x | (np.arange(N) << k)] = B[yp, :]
    return A % p
def tomat(A, p):
    n, m = A.shape; return nmod_mat(n, m, [int(v) for v in A.ravel()], p)
def ring_bounds(k, S1, S2):
    msk = 2**k - 1
    s1 = set(S1); s2 = set(S2)
    s12 = set(((t1 & msk) ^ (t2 >> k)) | (((t1 >> k) ^ (t2 & msk)) << k) for t1 in s1 for t2 in s2)
    a, b, c = len(s1), len(s2), len(s12)
    return min(a*b, a*c, b*c), math.sqrt(a*b*c), (a, b, c)
k = int(sys.argv[1]); N = 2**k; n = N*N
fr, beta = majorana_windows(k)
Mm = Mperm(k)
res = []
for i1, d1 in enumerate(fr):
    for i2, d2 in enumerate(fr):
        if len(sys.argv) > 2 and i1 != i2: continue
        rk = []
        for p in (P1, P2):
            A12 = tomat(build(k, d1['phiL'], 12, p), p); A23 = tomat(build(k, d2['phiL'], 23, p), p)
            s = (N*N*d1['L']*d2['L']) % p
            E = tomat((s*Mm) % p, p) - A12*A23
            rk.append(E.rank())
        up, lo, sizes = ring_bounds(k, d1['supp'], d2['supp'])
        r = rk[0]
        u = math.log(N*r)/(3*math.log(N)) if r > 0 else 0
        v_up = math.log(up)/(3*math.log(N)); v_lo = math.log(lo)/(3*math.log(N))
        print(f"k={k} W1={d1['win']}(q={d1['q']},p={d1['p']}) W2={d2['win']}(q={d2['q']},p={d2['p']}) rank(M-C')={rk} "
              f"[2p bound {d1['p']+d2['p']}] sigma-rank={N*r} ring in [{lo:.1f},{up}] proj={sizes} u={u:.3f} v in [{v_lo:.3f},{v_up:.3f}]", flush=True)
        res.append(dict(k=k, w1=d1['win'], w2=d2['win'], q1=d1['q'], q2=d2['q'], p1=d1['p'], p2=d2['p'], rank=rk, ring_up=up, ring_lo=lo, proj=sizes))
json.dump(res, open(f'composite_k{k}.json', 'w'))
