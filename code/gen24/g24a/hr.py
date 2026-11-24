# Hrubes beta-construction (TR24-026 Thm 9/12), exact, in our (p,q) setting.
import numpy as np, scipy.sparse as sp, itertools, math, random, sys, time
from fractions import Fraction
from math import comb
from cliff import family, rho
import ev
def alphas(k):
    """alpha_c (c = |C|) of the multilinear f with f(A) = sum_{C subset A} alpha_C; f(|a|=k)=1, f(|a|<k, =k mod 2)=0"""
    roots = [i for i in range(k) if (i - k) % 2 == 0]
    g = lambda w: math.prod(w - i for i in roots)
    c0 = Fraction(1, g(k))
    G = [c0 * g(w) for w in range(k + 1)]
    a = [sum((-1) ** (c - j) * comb(c, j) * G[j] for j in range(c + 1)) for c in range(k + 1)]
    d = k // 2
    assert all(x == 0 for x in a[d + 1:]), a
    return a[:d + 1]
def eA(E, A, n):
    M = np.eye(n, dtype=np.int64)
    for i in A: M = M @ E[i]
    return M
def setup(r, k, N=None, order='lex', seed=0):
    E = family(r); n = 2 ** r; t = len(E)
    sets = list(itertools.combinations(range(t), k))
    if order == 'rand': random.Random(seed).shuffle(sets)
    if N is None: N = n
    assert len(sets) >= N
    fam = sets[:N]
    al = alphas(k); d = k // 2
    den = math.lcm(*[a.denominator for a in al])
    Cs = [C for c in range(d + 1) for C in itertools.combinations(range(t), c) if al[c] != 0]
    Cs = [C for C in Cs if any(set(C) <= set(A) for A in fam)]       # pieces that actually occur
    eAs = [eA(E, A, n) for A in fam]
    return dict(E=E, n=n, t=t, k=k, d=d, N=N, fam=fam, al=al, den=den, Cs=Cs, eAs=eAs)
def H(S, a, tilde):
    """H_A (or H~_A) as n x (n*|Cs|) integer matrix (alpha scaled by den)"""
    n = S['n']; A = set(S['fam'][a]); cols = []
    for C in S['Cs']:
        if set(C) <= A:
            w = 1 if tilde else int(S['al'][len(C)] * S['den'])
            cols.append(w * S['eAs'][a])
        else: cols.append(np.zeros((n, n), dtype=np.int64))
    return np.hstack(cols)
def check_HR(S):
    N, n, den = S['N'], S['n'], S['den']
    Hs = [H(S, a, False) for a in range(N)]; Ht = [H(S, a, True) for a in range(N)]
    for a in range(N):
        if not (Hs[a] @ Ht[a].T == den * np.eye(n, dtype=np.int64)).all(): return False
        for b in range(a):
            if (Hs[a] @ Ht[b].T + Hs[b] @ Ht[a].T).any(): return False
    return True
def terms(S, pieces=None):
    """C-side terms (G, G') for ev: G(x=A,u=j) = e_A[j,r][C subset A], G'(y=i,v=B) = den*alpha_C e_B[i,r][C subset B];
    ev.build(N, terms, d=2*den) = 2den X - sum G (x) G'."""
    N, n = S['N'], S['n']
    if pieces is None: pieces = range(len(S['Cs']))
    out = []
    for ci in pieces:
        C = set(S['Cs'][ci]); w = int(S['al'][len(C)] * S['den'])
        mem = np.array([C <= set(A) for A in S['fam']], dtype=np.int64)
        for r_ in range(n):
            G = np.zeros((N, n), dtype=np.int64); Gp = np.zeros((n, N), dtype=np.int64)
            for a in range(N):
                if mem[a]:
                    G[a, :] = S['eAs'][a][:, r_]
                    Gp[:, a] = w * S['eAs'][a][:, r_]
            out.append((G, Gp))
    return out
def crank(S, tm):
    """exact rank of the C-flattening (x,u)|(y,v) of sum_l G_l (x) G'_l  = rank(U V^T), U = [vec G], V = [vec G']"""
    if not tm: return 0
    U = sp.csr_matrix(np.array([G.ravel() for G, _ in tm]).T); V = sp.csr_matrix(np.array([Gp.ravel() for _, Gp in tm]).T)
    return ev.rank(sp.csr_matrix(U @ V.T))[0]
def prank(S, tm):
    return ev.rank(ev.build(S['N'], tm, d=2 * S['den']))[0]
if __name__ == "__main__":
    r, k = int(sys.argv[1]), int(sys.argv[2])
    S = setup(r, k)
    t0 = time.time()
    print(f"n=N={S['n']} t={S['t']} k={k} d={S['d']} alphas={S['al']} pieces={len(S['Cs'])} (C(t,<=d)={sum(comb(S['t'],c) for c in range(S['d']+1))})")
    print("HR conditions exact:", check_HR(S), f"{time.time()-t0:.1f}s")
    tm = terms(S)
    print("full: q (#terms) =", len(tm), " q exact =", crank(S, tm), " p = rank(X - C) =", prank(S, tm), f"{time.time()-t0:.1f}s")
