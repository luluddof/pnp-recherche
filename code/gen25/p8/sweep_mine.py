# G25-P8 independent check of the sweep lemma (not derived from p2/sweep_check.py).
# (K,lam)-clique (paper Def. 3): one lam-bit string per edge {w,w'} of K_K; value = sum over all edge strings of
# prod_w W_w(strings on the K-1 edges at w, ordered by neighbour index).
import itertools, random, sys
import numpy as np

def edges_at(w, K):
    return [tuple(sorted((w, v))) for v in range(K) if v != w]

def brute(K, lam, W):
    E = list(itertools.combinations(range(K), 2)); D = 1 << lam
    idx = {e: i for i, e in enumerate(E)}
    # all assignments as a grid (vectorised)
    grids = np.indices([D] * len(E)).reshape(len(E), -1)
    tot = None
    for w in range(K):
        cols = [grids[idx[e]] for e in edges_at(w, K)]
        flat = np.ravel_multi_index(cols, [D] * (K - 1))
        vals = W[w].reshape(-1)[flat].astype(object)
        tot = vals if tot is None else tot * vals
    return int(sum(tot))

def sweep(K, lam, W):
    """eliminate nodes 0..K-1 in order; state: dict over assignments of the cut edges {u<=j<v}."""
    D = 1 << lam
    cut = []                    # current cut edges, list of (u,v)
    F = {(): 1}
    max_enum = 0; terms = 0
    for j in range(K):
        Dj = [e for e in cut if e[1] == j]          # edges into node j from earlier nodes (summed out)
        rest = [e for e in cut if e[1] != j]
        Nj = [(j, v) for v in range(j + 1, K)]      # new edges from j to later nodes
        newcut = rest + Nj
        max_enum = max(max_enum, len(cut) + len(Nj))
        G = {}
        for sig, val in F.items():
            a = dict(zip(cut, sig))
            for nv in itertools.product(range(D), repeat=len(Nj)):
                b = dict(a); b.update(zip(Nj, nv))
                arg = tuple(b[e] for e in edges_at(j, K))
                key = tuple(b[e] for e in newcut)
                G[key] = G.get(key, 0) + val * int(W[j][arg])
                terms += 1
        F = G; cut = newcut
    assert cut == [] and list(F.keys()) == [()]
    return F[()], max_enum, terms

def s_formula(K): return (K + 1) ** 2 // 4 - 1

if __name__ == "__main__":
    rng = np.random.default_rng(12345)
    ok = True
    for K, lam, reps in [(2, 1, 3), (2, 3, 3), (3, 1, 5), (3, 2, 5), (3, 3, 3), (4, 1, 5), (4, 2, 4), (5, 1, 5), (5, 2, 2), (6, 1, 3), (7, 1, 1)]:
        for r in range(reps):
            D = 1 << lam
            W = [rng.integers(-50, 51, size=[D] * (K - 1)) for _ in range(K)]
            b = brute(K, lam, W); s, me, tm = sweep(K, lam, W)
            good = (b == s) and (me == s_formula(K))
            ok &= good
            print(f"K={K} lam={lam} brute={b} sweep={s} eq={b==s} max_enum_edges={me} s(K)={s_formula(K)} "
                  f"floor(K^2/4)={K*K//4} K(K-1)/2={K*(K-1)//2} sweep_terms={tm} bound K*2^(lam s)={K*2**(lam*s_formula(K))}"
                  f" brute_terms={2**(lam*K*(K-1)//2)}")
            ok &= tm <= K * 2 ** (lam * s_formula(K))
    print("ALL OK" if ok else "FAIL")
