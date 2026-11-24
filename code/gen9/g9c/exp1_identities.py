"""G9c exp1: sanity checks of the class-restricted cut-off (Thm A: no checks; Thm B: checks + mod 2^nu).
Brute force over all subsets of a small graph.
(a) Thm A: for each size t, with nbar = max_{|S|=t} e(S), sum_{|S|=t} sum_{j<=nbar} (-1)^j C(e(S),j) = i_t(G);
    with nbar-1 the identity must (generically) fail.
(b) Thm B: class C_t = {|S|=t}, small cut-off nbar, random sparse checks until
    X = sum_{S in C_t} C(e(S), nbar+1) prod_j [|V_j(S)| even] = 0; then
    F = sum_{S in C_t} eta(S) * sum_{(Pi_j): sum|Pi_j| <= M'} prod_j sigma_{|Pi_j|}   (mod 2^nu)
    with nu = ceil(log2 C(n,t)) + 1, M' = nu - 1 + r  must equal i_t(G).  (Lemma 12 of the paper on a class.)
"""
import itertools, math, random, sys
import networkx as nx
from math import comb

def run(n, d, seed):
    rng = random.Random(seed)
    G = nx.random_regular_graph(d, n, seed=seed)
    E = list(G.edges()); m = len(E)
    emask = [(1 << u) | (1 << v) for u, v in E]
    # e(S) and violated-edge bitmask for all S
    viol = [0] * (1 << n)
    for S in range(1 << n):
        vm = 0
        for j, em in enumerate(emask):
            if S & em == em:
                vm |= 1 << j
        viol[S] = vm
    bysize = [[] for _ in range(n + 1)]
    for S in range(1 << n):
        bysize[S.bit_count()].append(S)
    it = [sum(1 for S in bysize[t] if viol[S] == 0) for t in range(n + 1)]
    alpha = max(t for t in range(n + 1) if it[t] > 0)
    out = dict(n=n, d=d, m=m, alpha=alpha, i=sum(it))
    # (a) Theorem A
    okA, failA, rowsA = True, 0, []
    for t in range(0, alpha + 2):
        es = [viol[S].bit_count() for S in bysize[t]]
        nb = max(es)
        def trunc(nbar):
            return sum(sum((-1) ** j * comb(e, j) for j in range(0, min(e, nbar) + 1)) for e in es)
        okA &= (trunc(nb) == it[t])
        if nb >= 1 and trunc(nb - 1) != it[t]:
            failA += 1
        costA = 0 if nb == 0 else math.log2(comb(m, nb))
        rowsA.append((t, it[t], nb, round(costA, 1)))
    out['thmA_ok'] = okA; out['thmA_fail_with_nbar_minus_1'] = failA; out['rowsA(t,i_t,nbar,logC(m,nbar))'] = rowsA
    # (b) Theorem B on classes t (use nbar = 1 or 2)
    sig = lambda i: 1 if i == 0 else (-1) ** i * 2 ** (i - 1)
    rowsB = []
    okB = True
    for t in range(2, alpha + 2):
        for nbar in (1, 2):
            cls = bysize[t]
            h = math.log2(comb(n, t)); nu = math.ceil(h) + 1
            checks = []
            alive = [S for S in cls if viol[S].bit_count() > nbar]
            X = sum(comb(viol[S].bit_count(), nbar + 1) for S in alive)
            X0 = X
            while X > 0:
                best = None
                for _ in range(40):
                    rate = min(0.5, 2.0 / max(nbar, 1) / 2)
                    c = 0
                    for j in range(m):
                        if rng.random() < rate:
                            c |= 1 << j
                    Xc = sum(comb(viol[S].bit_count(), nbar + 1) for S in alive if (viol[S] & c).bit_count() % 2 == 0)
                    if best is None or Xc < best[0]:
                        best = (Xc, c)
                X, c = best
                checks.append(c)
                alive = [S for S in alive if (viol[S] & c).bit_count() % 2 == 0]
            r = len(checks)
            Mp = nu - 1 + r
            F = 0; maxmarks = 0
            for S in cls:
                e = viol[S].bit_count()
                eta = sum((-1) ** j * comb(e, j) for j in range(0, min(e, nbar) + 1))
                if eta == 0:
                    continue
                # truncated product over checks: polynomial in z up to degree Mp
                poly = [1] + [0] * Mp
                for c in checks:
                    q = (viol[S] & c).bit_count()
                    new = [0] * (Mp + 1)
                    for a, ca in enumerate(poly):
                        if ca == 0: continue
                        for i in range(0, min(q, Mp - a) + 1):
                            new[a + i] += ca * sig(i) * comb(q, i)
                    poly = new
                F += eta * sum(poly)
            ok = (F % (1 << nu)) == it[t]
            okB &= ok
            rowsB.append(dict(t=t, nbar=nbar, h=round(h, 1), nu=nu, X0=X0, r=r, Mprime=Mp, F_mod=F % (1 << nu), i_t=it[t], ok=ok))
    out['thmB_ok'] = okB
    out['rowsB'] = rowsB
    return out

if __name__ == '__main__':
    for (n, d, seed) in [(12, 3, 1), (14, 4, 2), (16, 4, 3), (16, 6, 4)]:
        o = run(n, d, seed)
        print({k: v for k, v in o.items() if k not in ('rowsB',)})
        for rrow in o['rowsB']:
            print('   ', rrow)
        sys.stdout.flush()
