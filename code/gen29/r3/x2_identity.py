# G29-R3 check of X2's composite identity: sigma = S12 S23, SWAP = R + K  =>  sigma - K12 K23 = R12 S23 + K12 R23,
# rank <= 2pN, K12 K23 = ring/MPO with bonds (q,q,1) (flattening ranks site1|23 <= q, site3|12 <= q, site2|13 <= q^2).
# Own code: (a) random integer K of osr q (exact over Z), (b) Majorana (8,8)@8 and (6,19)@8 decompositions (mod primes).
import numpy as np, flint, itertools
rng = np.random.default_rng(5)
def rkZ(M): return flint.fmpz_mat([[int(v) for v in r] for r in M]).rank()
def rkP(M, p): return flint.nmod_mat(M.shape[0], M.shape[1], [int(x) for x in (M % p).ravel()], p).rank()
def swap(N):
    S = np.zeros((N*N, N*N), dtype=np.int64)
    for i in range(N):
        for j in range(N): S[j*N+i, i*N+j] = 1   # |i,j> -> |j,i>
    return S
def osr(K, N, rk): return rk(K.reshape(N, N, N, N).transpose(0, 2, 1, 3).reshape(N*N, N*N))
def flat(T, N, site, rk):  # T is N^3 x N^3 operator, rows (b1,b2,b3), cols (a1,a2,a3)
    T6 = T.reshape((N,)*6); other = [s for s in range(3) if s != site]
    return rk(T6.transpose(site, site+3, other[0], other[1], other[0]+3, other[1]+3).reshape(N*N, N**4))
def check(N, K, rk, mod=None, label=''):
    I = np.eye(N, dtype=np.int64); S = swap(N); R = S - K
    red = (lambda M: M % mod) if mod else (lambda M: M)
    p = rk(red(R)); q = osr(red(K), N, rk)
    S12, S23 = np.kron(S, I), np.kron(I, S); K12, K23 = np.kron(K, I), np.kron(I, K); R12, R23 = np.kron(R, I), np.kron(I, R)
    sig = S12 @ S23
    # sigma is the cyclic shift: which one?
    perm = {}
    for a in itertools.product(range(N), repeat=3):
        col = (a[0]*N+a[1])*N+a[2]; row = int(np.nonzero(sig[:, col])[0][0]); perm[a] = (row//(N*N), (row//N) % N, row % N)
    lhs = red(sig - K12 @ K23); rhs = red(R12 @ S23 + K12 @ R23)
    C = red(K12 @ K23)
    print(f'{label} N={N} (p,q)=({p},{q}): identity {np.array_equal(lhs, rhs)}; rank(sigma-K12K23)={rk(lhs)} <= 2pN={2*p*N}; '
          f'flattenings site1|23={flat(C,N,0,rk)}, site2|13={flat(C,N,1,rk)}, site3|12={flat(C,N,2,rk)} (bonds q,q,1 => <= {q},{q*q},{q}); '
          f'sigma maps (a1,a2,a3)->{perm[(0,1,2)] if N>=3 else perm[(0,1,1)]}', flush=True)
for N in (2, 3, 4):
    for q in (1, 2, 3):
        E = rng.integers(-2, 3, (q, N, N)); F = rng.integers(-2, 3, (q, N, N))
        K = sum(np.kron(E[l], F[l]) for l in range(q))
        check(N, K, rkZ, label='random')
# Majorana (8,8)@8 via Bell model, mod prime
P = 1000003; k = 3; N = 8; n2 = 64; msk = 7
POP = np.array([bin(i).count('1') for i in range(256)])
Q = lambda x: POP[(x & msk) & (x >> k)] & 1
B = lambda x, v: (POP[(x & msk) & (v >> k)] + POP[(x >> k) & (v & msk)]) & 1
def pauli(s):
    c, d = s & msk, s >> k; M = np.zeros((N, N), dtype=np.int64)
    for x in range(N): M[x ^ c, x] = -1 if bin(d & x).count('1') % 2 else 1
    return M
W = np.arange(n2)
for name, S_, kap in [('(6,19)', [0, 1, 8, 11, 25, 27], None)]:
    K = np.zeros((n2, n2), dtype=np.int64); inv2 = pow(2, -1, P)
    for s in S_: K = (K + (int(1-2*Q(s))*inv2 % P)*np.kron(pauli(s), pauli(s))) % P
    check(N, K, lambda M: rkP(M, P), mod=P, label='Bell'+name)
# K = I (osr 1, p = N(N-1)/2) for N = 2..5, exact over Z
for NN in (2, 3, 4, 5):
    check(NN, np.eye(NN*NN, dtype=np.int64), rkZ, label='K=I')
# Majorana (8,8)@8: window u in [1,2] of a JW-free frame, coefficients c(s) = N^-2 sum_v phi(v)(-1)^{B(s,v)} mod P
from fractions import Fraction
fr = None
for cand in itertools.combinations(range(1, n2), 6):
    if all(B(a, b) == 1 for a, b in itertools.combinations(cand, 2)): fr = list(cand); break
t = 0
for v in fr: t ^= v
fr.append(t)
Ja = sum(((Q(s) + B(W, s)) & 1) for s in fr); par = int(Ja[0] % 2); u = (Ja-par)//2; beta = 1-2*Q(W)
F = {j: int(beta[Ja % 4 == j][0]) for j in (par, par+2)}
Fu = lambda x: F[(2*x+par) % 4]
def f(x, win=(1, 2)):
    s = Fraction(0)
    for i in win:
        tt = Fraction(Fu(i))
        for j in win:
            if j != i: tt *= Fraction(x-j, i-j)
        s += tt
    return s
phi = [f(int(uu)) for uu in u]
K = np.zeros((n2, n2), dtype=np.int64)
for s in range(n2):
    c = sum((phi[v]*(1-2*int(B(s, v))) for v in range(n2)), Fraction(0))/n2
    if c: K = (K + (c.numerator % P)*pow(c.denominator % P, -1, P)*np.kron(pauli(s), pauli(s))) % P
check(N, K, lambda M: rkP(M, P), mod=P, label='Majorana window[1,2]')
