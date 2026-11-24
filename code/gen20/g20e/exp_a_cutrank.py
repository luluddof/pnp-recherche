"""Exp A: actual rank of the tensor carried across a variable bisection of a random 3-CNF.
For a bipartition (A,B) of the variables:
  crossing clauses c = clauses with variables on both sides; dA = A-variables in crossing clauses (similarly dB);
  W_A(alpha) = #assignments of A satisfying all clauses inside A, given alpha on dA (internal vars summed out);
  C(alpha,beta) = product of crossing clauses.
  The cut matrix  Mcut = diag(W_A) C diag(W_B);  Z = 1^T Mcut 1.
Quantities: bits = min(|dA|,|dB|) (what cut-conditioning / Claude's edge bits pay),  c (clause-rank bound 2^c),
  rank C (full), rank Mcut (= C restricted to reachable states), support sizes.
Ranks: exact over GF(p), p = 2^31-1 (lower bound on rank over Q; equal for 0/1 data w.h.p.), cross-checked with float SVD.
"""
import numpy as np, random, sys, itertools, networkx as nx

P = 2147483647

def rank_modp(M):
    M = (np.array(M, dtype=np.int64) % P).copy()
    r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i, c]:
                piv = i; break
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), P - 2, P)
        M[r] = (M[r] * inv) % P
        nz = np.nonzero(M[:, c])[0]
        for i in nz:
            if i != r:
                M[i] = (M[i] - M[i, c] * M[r]) % P
        r += 1
        if r == rows: break
    return r

def rank_float(M):
    if M.size == 0: return 0
    return int(np.linalg.matrix_rank(M.astype(float)))

def rand_cnf(n, m, rng):
    return [tuple((v, rng.random() < 0.5) for v in rng.sample(range(n), 3)) for _ in range(m)]

def sat(cl, asg):  # asg dict var->bool ; literal (v, neg)
    return any(asg[v] != neg for v, neg in cl)

def analyse(n, cls, A, B):
    A = sorted(A); B = sorted(B); Aset = set(A); Bset = set(B)
    inA = [c for c in cls if all(v in Aset for v, _ in c)]
    inB = [c for c in cls if all(v in Bset for v, _ in c)]
    cr = [c for c in cls if c not in inA and c not in inB]
    dA = sorted({v for c in cr for v, _ in c if v in Aset}); dB = sorted({v for c in cr for v, _ in c if v in Bset})
    def bitsarr(vars_, k):
        x = np.arange(2 ** k)[:, None] >> np.arange(k)[None, :] & 1
        return {v: x[:, i].astype(bool) for i, v in enumerate(vars_)}
    def W(side, inside, d):
        X = bitsarr(side, len(side))
        ok = np.ones(2 ** len(side), dtype=bool)
        for c in inside:
            ok &= np.any([X[v] != neg for v, neg in c], axis=0)
        a = np.zeros(2 ** len(side), dtype=np.int64)
        for i, v in enumerate(d): a += X[v].astype(np.int64) << i
        return np.bincount(a[ok], minlength=2 ** len(d)).astype(np.int64)
    WA = W(A, inA, dA); WB = W(B, inB, dB)
    XA = bitsarr(dA, len(dA)); XB = bitsarr(dB, len(dB))
    C = np.ones((2 ** len(dA), 2 ** len(dB)), dtype=np.int64)
    for c in cr:
        fa = np.ones(2 ** len(dA), dtype=bool); fb = np.ones(2 ** len(dB), dtype=bool)
        for v, neg in c:
            if v in Aset: fa &= (XA[v] == neg)
            else: fb &= (XB[v] == neg)
        C *= 1 - np.outer(fa, fb).astype(np.int64)
    sA = WA > 0; sB = WB > 0
    Mres = C[np.ix_(sA, sB)]
    Mcut = (WA[:, None] * C * WB[None, :])
    Z = int(Mcut.sum())
    rC = rank_modp(C) if max(C.shape) <= 256 else rank_float(C)
    rR = (rank_modp(Mres) if max(Mres.shape) <= 256 else rank_float(Mres)) if Mres.size else 0
    return dict(n=n, m=len(cls), cross=len(cr), dA=len(dA), dB=len(dB), suppA=int(sA.sum()), suppB=int(sB.sum()),
                rankC=rC, rankC_f=rank_float(C), rank_restr=rR, rank_restr_f=rank_float(Mres), Z=Z,
                bits=min(len(dA), len(dB)), log2rC=float(np.log2(max(rC, 1))), log2rR=float(np.log2(max(rR, 1))))

def bisect(n, cls, rng, how):
    G = nx.Graph(); G.add_nodes_from(range(n))
    for c in cls:
        vs = [v for v, _ in c]
        for i in range(3):
            for j in range(i + 1, 3):
                G.add_edge(vs[i], vs[j])
    if how == 'random':
        p = list(range(n)); rng.shuffle(p); return p[:n // 2], p[n // 2:]
    best = None
    for t in range(40):
        a, b = nx.algorithms.community.kernighan_lin_bisection(G, seed=rng.randrange(10 ** 9))
        # objective: boundary bits (what we must enumerate), not edges
        A, B = sorted(a), sorted(b)
        Aset = set(A)
        cr = [c for c in cls if not all(v in Aset for v, _ in c) and not all(v not in Aset for v, _ in c)]
        dA = {v for c in cr for v, _ in c if v in Aset}; dB = {v for c in cr for v, _ in c if v not in Aset}
        key = max(len(dA), len(dB))
        if best is None or key < best[0]: best = (key, A, B)
    return best[1], best[2]

if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    cfgs = [(20, 2.0, 'kl'), (20, 1.5, 'kl'), (24, 1.5, 'kl'), (28, 1.0, 'kl'), (32, 1.0, 'kl'), (20, 3.0, 'kl')]
    print('n m  how  cross |dA| |dB| suppA/2^dA suppB/2^dB  rankC (2^c bound)  rank_restr  log2rC log2rR bits  Z')
    for n, dens, how in cfgs:
        for rep in range(6):
            m = int(round(dens * n))
            cls = rand_cnf(n, m, rng)
            A, B = bisect(n, cls, rng, how)
            if max(len(A), len(B)) > 24: continue
            # skip too big boundaries
            Aset = set(A)
            cr = [c for c in cls if not all(v in Aset for v, _ in c) and not all(v not in Aset for v, _ in c)]
            dA = {v for c in cr for v, _ in c if v in Aset}; dB = {v for c in cr for v, _ in c if v not in Aset}
            if len(dA) > 10 or len(dB) > 10:
                print(n, m, how, 'boundary too large', len(dA), len(dB)); continue
            r = analyse(n, cls, A, B)
            print(r['n'], r['m'], how, r['cross'], r['dA'], r['dB'], f"{r['suppA']}/{2**r['dA']}", f"{r['suppB']}/{2**r['dB']}",
                  r['rankC'], f"(2^{r['cross']})", r['rank_restr'], f"{r['log2rC']:.2f}", f"{r['log2rR']:.2f}", r['bits'], r['Z'],
                  'float-check', r['rankC_f'] == r['rankC'], r['rank_restr_f'] == r['rank_restr'], flush=True)
