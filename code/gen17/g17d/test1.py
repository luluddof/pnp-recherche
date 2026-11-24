# G18c: independent exact tests for g14a sec.3 / sec.4(ii).  Arithmetic mod prime P (exact).
import itertools, random, sys
import numpy as np
from fractions import Fraction
from math import comb

def rref(M, P):
    M = M.copy() % P
    r = 0; piv = []
    rows, cols = M.shape
    for c in range(cols):
        if r == rows: break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        i = r + nz[0]
        if i != r: M[[r, i]] = M[[i, r]]
        M[r] = (M[r] * pow(int(M[r, c]), P - 2, P)) % P
        f = M[:, c].copy(); f[r] = 0
        idx = np.nonzero(f)[0]
        if len(idx):
            M[idx] = (M[idx] - np.outer(f[idx], M[r])) % P
        piv.append(c); r += 1
    return M[:r], piv

def rank(M, P):
    M = np.array(M, dtype=np.int64)
    if M.size == 0: return 0
    return len(rref(M, P)[1])

def nullspace(M, P):
    R, piv = rref(M, P)
    cols = M.shape[1]
    free = [c for c in range(cols) if c not in set(piv)]
    N = np.zeros((len(free), cols), dtype=np.int64)
    for i, f in enumerate(free):
        N[i, f] = 1
        for r, p in enumerate(piv):
            N[i, p] = (-R[r, f]) % P
    return N

def binp(l, m):  # Pr[Bin(l,1/2) <= m]
    if m < 0: return Fraction(0)
    return Fraction(sum(comb(l, j) for j in range(0, min(m, l) + 1)), 2 ** l)

def all_S(q, l):
    subs = list(range(2 ** q))
    return itertools.product(subs, repeat=l)

def wrow(lab, S):  # #{k: lab_k notin S_k}
    return sum(1 for k, s in enumerate(S) if not (s >> lab[k]) & 1)
def wcol(lab, S):  # #{k: lab_k in S_k}
    return sum(1 for k, s in enumerate(S) if (s >> lab[k]) & 1)

def build_instance(n, q, l, alpha, beta, degs, ranks, P, rng, taus, structuredU=False):
    """Random element of the solution space: A(t) with deg A[u,v] <= degs[u][v],
    A(tau_k) - U_k V_k block-diagonal by k-th symbol.  Returns coefficient list and E_k."""
    # unknown index
    idx = {}; cnt = 0
    for u in range(n):
        for v in range(n):
            for j in range(degs[u][v] + 1):
                idx[(u, v, j)] = cnt; cnt += 1
    nA = cnt
    vidx = {}
    for k in range(l):
        for s in range(ranks[k]):
            for v in range(n):
                vidx[(k, s, v)] = cnt; cnt += 1
    U = []
    for k in range(l):
        Uk = np.array([[rng.randrange(P) for _ in range(ranks[k])] for _ in range(n)], dtype=np.int64).reshape(n, ranks[k])
        if structuredU:  # each column supported on one k-th-symbol class of rows
            for s in range(ranks[k]):
                a = rng.randrange(q)
                for u in range(n):
                    if alpha[u][k] != a: Uk[u, s] = 0
        U.append(Uk)
    eqs = []
    for k in range(l):
        for u in range(n):
            for v in range(n):
                if alpha[u][k] != beta[v][k]:
                    row = np.zeros(cnt, dtype=np.int64)
                    for j in range(degs[u][v] + 1):
                        row[idx[(u, v, j)]] = pow(taus[k], j, P)
                    for s in range(ranks[k]):
                        row[vidx[(k, s, v)]] = (-U[k][u, s]) % P
                    eqs.append(row)
    M = np.array(eqs, dtype=np.int64) if eqs else np.zeros((0, cnt), dtype=np.int64)
    N = nullspace(M, P)
    if len(N) == 0: return None
    coef = np.zeros(cnt, dtype=np.int64)
    for i in range(len(N)):
        coef = (coef + rng.randrange(P) * N[i]) % P
    A0 = np.zeros((n, n), dtype=np.int64)
    for u in range(n):
        for v in range(n):
            A0[u, v] = coef[idx[(u, v, 0)]]
    E = []
    for k in range(l):
        V = np.zeros((ranks[k], n), dtype=np.int64)
        for s in range(ranks[k]):
            for v in range(n):
                V[s, v] = coef[vidx[(k, s, v)]]
        Ek = (U[k] @ V) % P
        # sanity: A(tau_k) - E_k is block diagonal
        for u in range(n):
            for v in range(n):
                val = sum(int(coef[idx[(u, v, j)]]) * pow(taus[k], j, P) for j in range(degs[u][v] + 1)) % P
                if alpha[u][k] != beta[v][k]:
                    assert (val - Ek[u, v]) % P == 0
        E.append(Ek)
    return A0, E, len(N)

def check(n, q, l, alpha, beta, pvec, qvec, d, ranks, P, rng, structuredU=False):
    """entry degree deg A[u,v] <= d + p_u + q_v.  Budget: w(u)+p_u <= M1, w'(v)+q_v <= M2, M1+M2+d = l-1."""
    taus = rng.sample(range(1, P), l)
    degs = [[d + pvec[u] + qvec[v] for v in range(n)] for u in range(n)]
    inst = build_instance(n, q, l, alpha, beta, degs, ranks, P, rng, taus, structuredU)
    if inst is None: return None
    A0, E, nd = inst
    tot = sum(rank(Ek, P) for Ek in E)
    L = l - 1 - d
    worst = 0; worst_ref = 0
    for S in all_S(q, l):
        wr = [wrow(alpha[u], S) + pvec[u] for u in range(n)]
        wc = [wcol(beta[v], S) + qvec[v] for v in range(n)]
        ref = 0
        for i in range(0, L + 1):
            R = [u for u in range(n) if wr[u] <= i]
            C = [v for v in range(n) if wc[v] <= L - i]
            rk = rank(A0[np.ix_(R, C)], P) if R and C else 0
            worst = max(worst, rk)
            ref += rk
        worst_ref = max(worst_ref, ref)
    return tot, worst, worst_ref, rank(A0, P), nd

if __name__ == "__main__":
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    # ---- Test A: averaging identity  E_S |G_S| = |dom pi| p^2
    okA = True
    for q in (2, 3):
        for l in (1, 2, 3, 4):
            if q ** l > 81 and False: continue
            for d in range(0, l):
                m = (l - d - 1) // 2
                n = 6
                alpha = [tuple(rng.randrange(q) for _ in range(l)) for _ in range(n)]
                beta = [tuple((a + 1 + rng.randrange(q - 1)) % q for a in alpha[u]) for u in range(n)]  # pi = id, differs everywhere
                totG = 0; cntS = 0
                for S in all_S(q, l):
                    totG += sum(1 for u in range(n) if wrow(alpha[u], S) <= m and wcol(beta[u], S) <= m)
                    cntS += 1
                if Fraction(totG, cntS) != n * binp(l, m) ** 2:
                    okA = False; print("AVERAGING FAILS", q, l, d)
    print("Test A (averaging identity, q=2,3, l<=4):", "OK" if okA else "FAIL")
    # ---- Test B: deterministic inequality, random + structured low-rank constructions
    viol = 0; runs = 0; tight = 0; best = []
    for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 60):
        q = rng.choice((2, 3)); l = rng.choice((2, 3, 4)); P = rng.choice((5, 7, 101))
        n = rng.choice((4, 6, 8)) if l < 4 else rng.choice((4, 6))
        d = rng.randrange(0, l)
        inj = rng.random() < 0.5
        alpha = [tuple(rng.randrange(q) for _ in range(l)) for _ in range(n)]
        mode = rng.choice(("iii", "free"))
        if mode == "iii":
            beta = [tuple((a + 1 + rng.randrange(q - 1)) % q for a in alpha[u]) for u in range(n)]
        else:
            beta = [tuple(rng.randrange(q) for _ in range(l)) for _ in range(n)]
        nonuni = rng.random() < 0.4
        pvec = [rng.choice((0, 0, 1)) if nonuni else 0 for _ in range(n)]
        qvec = [rng.choice((0, 0, 1)) if nonuni else 0 for _ in range(n)]
        ranks = [rng.choice((0, 1, 1, 2)) for _ in range(l)]
        res = check(n, q, l, alpha, beta, pvec, qvec, d, ranks, P, rng, structuredU=rng.random() < 0.5)
        if res is None: continue
        tot, worst, worst_ref, rkA, nd = res
        runs += 1
        if worst > tot or worst_ref > tot:
            viol += 1
            print("VIOLATION", dict(q=q, l=l, P=P, n=n, d=d, alpha=alpha, beta=beta, p=pvec, qv=qvec, ranks=ranks), res)
        if worst_ref == tot and tot > 0: tight += 1
        best.append((q, l, d, n, P, tot, worst, worst_ref, rkA))
    print("Test B: runs", runs, "violations", viol, "refined-bound tight in", tight)
    for b in best[:12]: print("  (q,l,d,n,P, sum rank E, max_S rank A0[R,C], max_S refined, rank A0) =", b)
