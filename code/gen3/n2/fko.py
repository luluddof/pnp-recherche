"""FKO-style even-cover certificates for random 3-CNF, measured against the spectral bound.
Identity: for x in {±1}^n, y_i = s_i x_{v_i} (true=+1), SAT x => sum_C (e1 - e2 + e3) = m.
 E1-E2 <= B12 := L + n*lmax(A)  (agent1 spectral, sound)
 E3 = sum_C tau_C, tau_C = s1 s2 s3 x_a x_b x_c ; V := #{C: tau_C=-1} => E3 = m-2V
 SAT => V <= B12/2.   Even cover K (every var even times): prod_{C in K} tau_C = prod_{C in K} s(C) =: pi_K.
 pi_K = -1 ("bad cover") => K contains an odd (>=1) number of C with tau_C=-1.
 Certificate: LP  min sum v_C  s.t. sum_{C in K} v_C >= 1 (bad K), 0<=v<=1  ; if LP* > B12/2 => UNSAT.
"""
import sys, time, math, random, itertools
from collections import deque, defaultdict
import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from scipy.optimize import linprog

def gen(n, m, rng):
    V = np.array([rng.sample(range(n), 3) for _ in range(m)])
    S = np.array([[rng.choice((1, -1)) for _ in range(3)] for _ in range(m)])
    return V, S

def spectral(n, V, S):
    m = len(V)
    c = np.zeros(n); A = np.zeros((n, n))
    for (a, b, d), s in zip(V, S):
        v = (a, b, d)
        for i in range(3): c[v[i]] += s[i]
        for i, j in ((0, 1), (0, 2), (1, 2)):
            A[v[i], v[j]] -= s[i]*s[j]/2; A[v[j], v[i]] -= s[i]*s[j]/2
    L = np.abs(c).sum(); Qb = n*max(0., np.linalg.eigvalsh(A)[-1])
    rows, cols, vals = [], [], []
    for v, s in zip(V, S):
        p = s[0]*s[1]*s[2]
        for i, j, k in itertools.permutations(range(3)):
            rows.append(v[i]); cols.append(v[j]*n + v[k]); vals.append(p/6)
    W = sp.csr_matrix((vals, (rows, cols)), shape=(n, n*n)); G = (W.T @ W).tocsr()
    trG = G.diagonal().sum(); Goff = G - sp.diags(G.diagonal())
    lg = eigsh(Goff.astype(float), k=1, which='LA', return_eigenvectors=False)[0]
    B3 = math.sqrt(n*(trG + n*n*max(0., lg)))
    return L + Qb, B3

def cycle_covers(n, V, roots, rng):
    """Even covers = cycles in H1: vertices = vars ∪ pairs, clause (a,b,c) gives edges a-{b,c}, b-{a,c}, c-{a,b}.
    Fundamental cycles of BFS trees (2-core only) from several roots. Poly time O(roots * |E|)."""
    m = len(V); pid = {}; edges = []
    for ci, (a, b, c) in enumerate(V):
        for u, p in ((a, (b, c)), (b, (a, c)), (c, (a, b))):
            key = (min(p), max(p))
            if key not in pid: pid[key] = n + len(pid)
            edges.append((u, pid[key], ci))
    N = n + len(pid); adj = defaultdict(list)
    for ei, (u, w, ci) in enumerate(edges):
        adj[u].append((w, ei)); adj[w].append((u, ei))
    deg = {v: len(adj[v]) for v in adj}; alive = set(adj); q = deque(v for v in adj if deg[v] < 2)
    while q:  # 2-core
        v = q.popleft()
        if v not in alive: continue
        alive.discard(v)
        for w, _ in adj[v]:
            if w in alive:
                deg[w] -= 1
                if deg[w] < 2: q.append(w)
    core = [v for v in alive if v < n]
    covers = set()
    if not core: return covers, len(alive)
    for r in rng.sample(core, min(roots, len(core))):
        par = {r: (None, None)}; dep = {r: 0}; dq = deque([r]); tree = set()
        while dq:
            v = dq.popleft()
            for w, ei in adj[v]:
                if w in alive and w not in par:
                    par[w] = (v, ei); dep[w] = dep[v]+1; tree.add(ei); dq.append(w)
        for ei, (u, w, ci) in enumerate(edges):
            if ei in tree or u not in par or w not in par: continue
            cl = defaultdict(int); cl[ci] ^= 1; a, b = u, w
            while a != b:
                if dep[a] < dep[b]: a, b = b, a
                pv, pe = par[a]; cl[edges[pe][2]] ^= 1; a = pv
            K = frozenset(k for k, t in cl.items() if t)
            if K and len(K) <= 40: covers.add(K)
    return covers, len(alive)

def check_even(K, V):
    cnt = defaultdict(int)
    for k in K:
        for v in V[k]: cnt[v] += 1
    return all(c % 2 == 0 for c in cnt.values())

def gauss_kernel(n, V, rng, tries=3):
    """GF(2) left kernel of m x n incidence; returns kernel dim and min weight among basis vectors over random row orders."""
    m = len(V); best = None; dim = None
    for _ in range(tries):
        order = list(range(m)); rng.shuffle(order)
        piv = {}  # pivot col -> (rowbits, combo bits)
        ws = []
        for i in order:
            r = 0
            for v in V[i]: r ^= 1 << int(v)
            comb = 1 << i
            while r:
                h = r.bit_length()-1
                if h in piv: pr, pc = piv[h]; r ^= pr; comb ^= pc
                else: piv[h] = (r, comb); break
            if r == 0: ws.append(bin(comb).count('1'))
        dim = len(ws); mn = min(ws) if ws else None
        if mn is not None and (best is None or mn < best): best = mn
        avg = sum(ws)/len(ws) if ws else None
    return dim, best, avg

def lp_bound(m, bad):
    if not bad: return 0.0
    rows, cols = [], []
    for j, K in enumerate(bad):
        for k in K: rows.append(j); cols.append(k)
    Aub = -sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(bad), m))
    res = linprog(np.ones(m), A_ub=Aub, b_ub=-np.ones(len(bad)), bounds=(0, 1), method='highs-ds')
    return res.fun

CAP = 12000
def run(n, m, seed, roots=30, do_gauss=True):
    rng = random.Random(seed); V, S = gen(n, m, rng)
    t = time.time(); B12, B3 = spectral(n, V, S); tsp = time.time()-t
    t = time.time(); covs, core = cycle_covers(n, V, roots, rng); tcy = time.time()-t
    sgn = np.prod(S, axis=1)
    covs = [K for K in covs if check_even(K, V)]
    bad = sorted([K for K in covs if np.prod([sgn[k] for k in K]) == -1], key=len)
    nbad_all = len(bad); bad = bad[:CAP]
    t = time.time(); lp = lp_bound(m, bad); tlp = time.time()-t
    sizes = sorted(len(K) for K in covs)
    g = gauss_kernel(n, V, rng, 1) if do_gauss else (None, None, None)
    return dict(n=n, m=m, a15=m/n**1.5, a14=m/n**1.4, B12=B12, B3=B3, spec_ok=B12+B3 < m,
                core=core, ncov=len(covs), nbad=nbad_all, minsz=sizes[0] if sizes else None,
                medsz=sizes[len(sizes)//2] if sizes else None, lp=lp, need=B12/2, cov_ok=lp > B12/2,
                kdim=g[0], gmin=g[1], gavg=g[2], tsp=tsp, tcy=tcy, tlp=tlp)

if __name__ == '__main__':
    n = int(sys.argv[1]); alphas = [float(a) for a in sys.argv[2].split(',')]
    for a in alphas:
        m = int(a*n**1.5)
        r = run(n, m, seed=1000*n+int(10*a), roots=(30 if n<=100 else 8), do_gauss=(n <= 200))
        print(' '.join(f'{k}={v:.3g}' if isinstance(v, float) else f'{k}={v}' for k, v in r.items()), flush=True)
