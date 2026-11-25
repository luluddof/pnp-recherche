# Independent check of the singlet-number family, built in the computational basis (no Bell/Pauli formulas used).
# Space C^N (x) C^N, N=2^k, index x*N+y.  SWAP_i swaps bit i of x with bit i of y.  L = sum_i SWAP_i, M = (k - L)/2.
# K = f(M), f = Lagrange interpolant of (-1)^m on window W.  p = rank(SWAP-K), q = osr(K) = rank of realignment.
import numpy as np, sys, time
from fractions import Fraction
from math import comb, lcm
PR = [1000000007, 998244353]
def rank_mod(A, p):
    A = A % p; n, m = A.shape; r = 0
    for c in range(m):
        if r == n: break
        nz = np.flatnonzero(A[r:, c])
        if nz.size == 0: continue
        i = r + nz[0]
        if i != r: A[[r, i]] = A[[i, r]]
        A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
        col = A[:, c].copy(); col[r] = 0
        idx = np.flatnonzero(col)
        if idx.size: A[idx] = (A[idx] - np.outer(col[idx], A[r]) % p) % p
        r += 1
    return r
def poly_in_L(k, W):
    # g(L) = f((k-L)/2) as polynomial in L with Fraction coeffs (via solving Vandermonde exactly on L-values)
    # f interpolates (-1)^m on W (deg |W|-1); L-value for m is k-2m.
    pts = [(k - 2*m, Fraction((-1)**m)) for m in W]
    d = len(pts) - 1
    # Newton -> monomial
    coeffs = [Fraction(0)]*(d+1)
    for i, (xi, yi) in enumerate(pts):
        basis = [Fraction(1)]; den = Fraction(1)
        for j, (xj, _) in enumerate(pts):
            if j == i: continue
            basis = [Fraction(0)] + basis
            for t in range(len(basis)-1): basis[t] -= xj * basis[t+1]
            den *= (xi - xj)
        for t in range(d+1): coeffs[t] += yi * basis[t] / den
    return coeffs
def build(k):
    N = 2**k; n = N*N
    x = np.repeat(np.arange(N), N); y = np.tile(np.arange(N), N)
    Ls = np.zeros((n, n), dtype=np.float64)
    for i in range(k):
        bx = (x >> i) & 1; by = (y >> i) & 1
        x2 = x ^ (bx << i) ^ (by << i); y2 = y ^ (by << i) ^ (bx << i)
        Ls[x2*N + y2, np.arange(n)] += 1
    SW = np.zeros((n, n), dtype=np.int64); SW[y*N + x, np.arange(n)] = 1
    return N, n, Ls, SW
def realign(K, N):
    # K[(x,y),(u,v)] -> R[(x,u),(y,v)]
    return K.reshape(N, N, N, N).transpose(0, 2, 1, 3).reshape(N*N, N*N)
def run(k, windows):
    N, n, Ls, SW = build(k)
    pows = [np.eye(n)]
    for j in range(k): pows.append(pows[-1] @ Ls)   # exact: entries <= k^j
    assert np.abs(pows[-1]).max() < 2**52
    pows = [P.astype(np.int64) for P in pows]
    # sanity: SWAP == prod(1-2E_i) == ... check (-1)^M via interpolation over full window gives SWAP exactly
    out = []
    for W in windows:
        c = poly_in_L(k, W); D = 1
        for t in c: D = lcm(D, t.denominator)
        ci = [int(t * D) for t in c]
        res = {}
        for p in PR:
            K = np.zeros((n, n), dtype=np.int64)
            for j, cj in enumerate(ci):
                if cj: K = (K + (cj % p) * (pows[j] % p)) % p
            A = (D % p) * SW % p - K
            r1 = rank_mod(A % p, p)
            r2 = rank_mod(realign(K, N) % p, p)
            res[p] = (r1, r2)
        mass = lambda m: comb(k, m) * 3**(k-m)
        p_f = sum(mass(m) for m in range(k+1) if m not in W)
        d = len(W) - 1
        q_f = sum(comb(k, u) * 3**u for u in range(d+1))
        a, b = res[PR[0]], res[PR[1]]
        ok = (a == b == (p_f, q_f))
        out.append((W, a, b, p_f, q_f, ok))
        print(f"k={k} N={N} W={W[0]}..{W[-1]} d={d}: rank(SWAP-K) mod p1,p2 = {a[0]},{b[0]}  osr = {a[1]},{b[1]}  "
              f"formula (p,q)=({p_f},{q_f}) {'OK' if ok else 'MISMATCH'}", flush=True)
    return out
if __name__ == "__main__":
    k = int(sys.argv[1]); mode = sys.argv[2] if len(sys.argv) > 2 else "all"
    wins = []
    for d in range(0, k+1):
        cand = [list(range(a, a+d+1)) for a in range(0, k-d+1)]
        if mode == "best":
            cand = [max(cand, key=lambda W: sum(comb(k, m)*3**(k-m) for m in W))]
        wins += cand
    t = time.time(); run(k, wins); print("time %.1fs" % (time.time() - t))
