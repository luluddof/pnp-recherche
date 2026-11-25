# Independent check of C2 Lemmas 2,3 (and improved threshold u in Lemma 3) for general d.
import itertools, sys, numpy as np
from math import factorial
rng = np.random.default_rng(7)
def partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield (); return
    for a in range(min(n, mx), 0, -1):
        for r in partitions(n-a, a): yield (a,) + r
def mn_char(shape, ctype):
    # Murnaghan-Nakayama via beta-sets
    if sum(shape) == 0: return 1
    r = ctype[0]; rest = ctype[1:]
    L = len(shape); beta = [shape[i] + L - 1 - i for i in range(L)]
    tot = 0
    for i, b in enumerate(beta):
        nb = b - r
        if nb < 0 or nb in beta: continue
        sign = (-1)**sum(1 for c in beta if nb < c < b)
        newb = sorted([c for c in beta if c != b] + [nb], reverse=True)
        Ln = len(newb); newshape = tuple(x for x in (newb[j] - (Ln-1-j) for j in range(Ln)) if x > 0)
        tot += sign * mn_char(newshape, rest)
    return tot
def cyc(pm):
    k = len(pm); seen = set(); L = []
    for i in range(k):
        if i in seen: continue
        j = i; c = 0
        while j not in seen: seen.add(j); j = pm[j]; c += 1
        L.append(c)
    return tuple(sorted(L, reverse=True))
def fdim(shape): return mn_char(shape, (1,)*sum(shape))
d, k = int(sys.argv[1]), int(sys.argv[2]); N = d**k
perms = list(itertools.permutations(range(k)))
idx = list(itertools.product(range(d), repeat=k)); pos = {t: i for i, t in enumerate(idx)}
def pmat(pm):
    M = np.zeros((N, N))
    for i, t in enumerate(idx):
        new = [None]*k
        for s in range(k): new[pm[s]] = t[s]
        M[pos[tuple(new)], i] = 1
    return M
Pm = {pm: pmat(pm) for pm in perms}
shapes = list(partitions(k))
ch = {(sh, pm): mn_char(sh, cyc(pm)) for sh in shapes for pm in perms}
Pi = {}
for sh in shapes:
    if len(sh) > d: continue
    Pi[sh] = sum(fdim(sh)/factorial(k)*ch[(sh, pm)]*Pm[pm] for pm in perms)
print("sector dims", {sh: round(np.trace(P)) for sh, P in Pi.items()}, "sum", round(sum(np.trace(P) for P in Pi.values())))
def iso(mu, A):  # conjugation-isotypic projection
    return sum(fdim(mu)/factorial(k)*ch[(mu, pm)]*(Pm[pm] @ A @ Pm[pm].T) for pm in perms)
def contains(a, b): return all((a[i] if i < len(a) else 0) >= x for i, x in enumerate(b))
def common_sub(nu, nup, size):
    return any(contains(nu, al) and contains(nup, al) for al in partitions(size))
# Lemma 2
bad = 0
for mu in shapes:
    if len(mu) > d*d: continue
    j = k - mu[0]
    A = iso(mu, rng.normal(size=(N, N)))
    if np.abs(A).max() < 1e-9: continue
    for nu in Pi:
        for nup in Pi:
            nz = np.abs(Pi[nu] @ A @ Pi[nup]).max() > 1e-9
            pred = common_sub(nu, nup, k-j)
            if nz and not pred: bad += 1; print("VIOLATION", mu, nu, nup)
            if nz and abs(nu[0]-nup[0]) > j: bad += 1; print("BAND VIOLATION", mu, nu, nup)
print("Lemma2 violations:", bad)
# Lemma 3 with threshold u (improved) and 2u
S = np.zeros((N*N, N*N))
for x in range(N):
    for y in range(N): S[y*N+x, x*N+y] = 1
def degop(u):
    return sum(iso(mu, rng.normal(size=(N, N))) for mu in shapes if mu[0] >= k-u and len(mu) <= d*d)
for u in range(0, k):
    nterms = 3
    K = sum(np.kron(degop(u), degop(u)) for _ in range(nterms))
    for thr in [u-1, u, 2*u]:
        if thr < 0: continue
        PY = sum(np.kron(Pi[a], Pi[b]) for a in Pi for b in Pi if a[0]-b[0] > thr)
        PX = sum(np.kron(Pi[a], Pi[b]) for a in Pi for b in Pi if b[0]-a[0] > thr)
        if isinstance(PY, int): print(f"u={u} thr={thr}: Y empty"); continue
        leak = np.abs(PX @ K @ PY).max(); dimY = round(np.trace(PY))
        r = np.linalg.matrix_rank(S - K, tol=1e-8)
        print(f"u={u} thr={thr}: |P_X K P_Y|max={leak:.1e} dimY={dimY} rank(SWAP-K)={r}")
