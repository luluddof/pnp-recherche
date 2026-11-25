# G29-R3: explicit-operator check of Majorana points. K = sum_s c(s) P_s (x) P_s built as an N^2 x N^2 matrix
# (P_(c,d) = X^c Z^d real Paulis), p := rank(SWAP - K), q := rank of the realigned K (= operator-Schmidt rank),
# both mod two primes.  c(s) = N^-2 sum_v phi(v) (-1)^{B(s,v)}, phi = window interpolant of F(u) on a RANDOM frame.
import sys, random
from fractions import Fraction
import numpy as np, flint
random.seed(int(sys.argv[2])); k = int(sys.argv[1]); N = 2**k; n2 = N*N; msk = N-1
POP = np.array([bin(i).count('1') for i in range(n2)], dtype=np.int64)
def Q(x): return POP[(x & msk) & (x >> k)] & 1
def B(x, v): return (POP[(x & msk) & (v >> k)] + POP[(x >> k) & (v & msk)]) & 1
W = np.arange(n2, dtype=np.int64)
def frame():
    for _ in range(1000):
        S = []
        cand = list(range(1, n2)); random.shuffle(cand)
        for v in cand:
            if all(B(np.int64(v), np.int64(s)) == 1 for s in S): S.append(v)
            if len(S) == 2*k: break
        if len(S) == 2*k:
            t = 0
            for v in S: t ^= v
            return S + [t]
fr = frame()
Ja = sum(((Q(np.int64(s)) + B(W, np.int64(s))) & 1) for s in fr)
beta = 1 - 2*Q(W); par = int(Ja[0] % 2); u = (Ja-par)//2; umax = (2*k+1-par)//2
F = {}
for j in (par, par+2): F[j] = int(beta[Ja % 4 == j][0]); assert len(set(beta[Ja % 4 == j])) == 1
Fu = lambda x: F[(2*x+par) % 4]
def pauli(s):
    c, d = s & msk, s >> k; M = np.zeros((N, N), dtype=np.int64)
    for x in range(N): M[x ^ c, x] = -1 if bin(d & x).count('1') % 2 else 1
    return M
PA = [pauli(s) for s in range(n2)]
SW = np.zeros((n2, n2), dtype=np.int64)
for i in range(N):
    for j in range(N): SW[i*N+j, j*N+i] = 1
Bmat = np.array([[B(np.int64(s), np.int64(v)) for v in range(n2)] for s in range(n2)]) if n2 <= 256 else None
def coeffs(phi):  # exact rational c(s)
    out = {}
    for s in range(n2):
        sg = 1 - 2*B(np.int64(s), W)
        tot = sum((Fraction(phi[v])*int(sg[v]) for v in range(n2)), Fraction(0)) / n2
        if tot != 0: out[s] = tot
    return out
def rank_mod(M, p):
    return flint.nmod_mat(M.shape[0], M.shape[1], [int(x) for x in (M % p).ravel()], p).rank()
primes = [1000003, 998244353]
wins = [tuple(map(int, w.split('-'))) for w in sys.argv[3].split(',')] if len(sys.argv) > 3 else [(lo, hi) for lo in range(umax+1) for hi in range(lo, umax+1)]
for lo, hi in wins:
    win = range(lo, hi+1)
    def f(x):
        s = Fraction(0)
        for i in win:
            t = Fraction(Fu(i))
            for j in win:
                if j != i: t *= Fraction(x-j, i-j)
            s += t
        return s
    fv = [f(x) for x in range(umax+1)]
    phi = [fv[uu] for uu in u]
    C = coeffs(phi)
    res = []
    for p in primes:
        K = np.zeros((n2, n2), dtype=np.int64)
        for s, cs in C.items():
            cm = cs.numerator % p * pow(cs.denominator % p, -1, p) % p
            K = (K + cm * np.kron(PA[s], PA[s])) % p
        rp = rank_mod((SW - K) % p, p)
        Rl = K.reshape(N, N, N, N).transpose(0, 2, 1, 3).reshape(n2, n2)
        rq = rank_mod(Rl, p)
        res.append((rq, rp))
    pw = sum(1 for uu in u if fv[uu] != Fu(uu))
    print(f'N={N} window u in [{lo},{hi}]: (q,p) explicit mod {primes}: {res} | #nonzero coeffs {len(C)} | p from phi {pw}', flush=True)
