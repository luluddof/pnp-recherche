# G-invariant crossing model: K = sum_O c_O 1_O over orbits O of G on X^4 (index (u,v,x,y), K[(u,v),(x,y)]).
# p = rank(SWAP-K), q = rank(R(K)), R(K)[(u,x),(v,y)] = K[(u,v),(x,y)].  Search block-vanishing patterns on isotypic
# components (central idempotents of A = End_G(V(x)V)).
import numpy as np, itertools, sys, time, json
from math import comb

def wreath_gens(ns):
    """generators of S_{n1} wr S_{n2} wr ... acting on leaves of tree with branching ns[0] (top) ... ns[-1] (bottom).
    Leaves indexed by digits (a_0,...,a_{L-1}), a_j in [ns[j]]; generator: at depth j, for each prefix, transpositions/cycle
    of digit j (acting only under one fixed prefix suffices, combined with top-level ones to generate)."""
    L = len(ns); leaves = list(itertools.product(*[range(n) for n in ns])); idx = {t: i for i, t in enumerate(leaves)}
    gens = []
    for j in range(L):
        n = ns[j]
        pre = tuple([0]*j)  # act under prefix 0...0 only
        for perm in ([(1, 0) + tuple(range(2, n))] + ([tuple(list(range(1, n)) + [0])] if n > 2 else [])):
            g = []
            for t in leaves:
                if t[:j] == pre:
                    t2 = t[:j] + (perm[t[j]],) + t[j+1:]
                else:
                    t2 = t
                g.append(idx[t2])
            gens.append(np.array(g))
    return gens, len(leaves)

def orbits4(gens, N):
    M = N**4
    parent = np.arange(M)
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    idx = np.arange(M); u = idx // N**3; v = (idx // N**2) % N; x = (idx // N) % N; y = idx % N
    for g in gens:
        img = ((g[u] * N + g[v]) * N + g[x]) * N + g[y]
        for a, b in zip(idx, img):
            ra, rb = find(a), find(b)
            if ra != rb: parent[ra] = rb
    roots = np.array([find(a) for a in idx])
    _, lab = np.unique(roots, return_inverse=True)
    return lab  # orbit label of (u,v,x,y)

def build(ns):
    gens, N = wreath_gens(ns)
    lab = orbits4(gens, N); r = lab.max() + 1
    D = N*N
    lab2 = lab.reshape(D, D)                     # rows (u,v), cols (x,y)
    # realignment permutation of orbits: orbit of (u,v,x,y) -> orbit of (u,x,v,y)
    L4 = lab.reshape(N, N, N, N)
    Rlab = L4.transpose(0, 2, 1, 3)              # Rlab[u,x,v,y] = lab[u,v,x,y]
    # orbit O as R-matrix: entries (u,x),(v,y) with lab(u,v,x,y)=O ; as an element of A it is orbit lab of (u,x,v,y) indexing
    piR = np.zeros(r, dtype=int)
    flat = lab.copy(); idx = np.arange(N**4)
    u = idx // N**3; v = (idx // N**2) % N; x = (idx // N) % N; y = idx % N
    other = ((u * N + x) * N + v) * N + y
    piR[lab] = lab[other]  # orbit of (u,v,x,y) maps to orbit of tuple (u,x,v,y) read as (u',v',x',y')
    reps = np.zeros(r, dtype=int); reps[lab[::-1]] = idx[::-1]  # some representative
    sizes = np.bincount(lab)
    return dict(N=N, r=r, lab2=lab2, piR=piR, reps=reps, sizes=sizes)

def mat(c, S):
    return c[S['lab2']]

def coeffs(Mx, S):
    """coefficients of an element of A in orbit basis (read at representatives)"""
    D = S['N']**2; f = Mx.reshape(-1)
    return f[S['reps']]

def center_idempotents(S, rng):
    r = S['r']; D = S['N']**2
    B = [mat(np.eye(r)[o], S) for o in range(r)]
    A1 = mat(rng.standard_normal(r), S); A2 = mat(rng.standard_normal(r), S)
    rows = []
    for A in (A1, A2):
        Mc = np.array([coeffs(Bo @ A - A @ Bo, S) for Bo in B]).T  # r x r
        rows.append(Mc)
    Mall = np.vstack(rows)
    u, s, vt = np.linalg.svd(Mall)
    nul = vt[np.sum(s > 1e-8 * s[0]):].T
    z = nul @ rng.standard_normal(nul.shape[1])
    Z = mat(z, S); Z = Z + Z.T
    ev, evec = np.linalg.eigh(Z)
    # cluster
    groups = []; cur = [0]
    for i in range(1, D):
        if abs(ev[i] - ev[cur[-1]]) < 1e-6 * max(1, abs(ev).max()): cur.append(i)
        else: groups.append(cur); cur = [i]
    groups.append(cur)
    E = [evec[:, g] @ evec[:, g].T for g in groups]
    return E, nul.shape[1]

def rk(M):
    s = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(s > 1e-7 * max(1.0, s[0])))

def search(ns, maxsub=12, seed=0, out=None):
    rng = np.random.default_rng(seed); t0 = time.time()
    S = build(ns); N, r = S['N'], S['r']; D = N*N
    E0, zc = center_idempotents(S, rng)
    SW = np.zeros((D, D)); a = np.arange(D); SW[(a % N) * N + a // N, a] = 1
    E = []
    for e in E0:
        for sg in (1, -1):
            f = e @ (np.eye(D) + sg * SW) / 2
            if np.trace(f) > 0.5: E.append(f)
    w = [int(round(np.trace(e))) for e in E]
    # linear maps c -> coeffs(K E_rho), c -> coeffs(R(K) E_rho)
    B = [mat(np.eye(r)[o], S) for o in range(r)]
    LR = [np.array([coeffs(Bo @ e, S) for Bo in B]).T for e in E]          # r x r
    sR = [coeffs(SW @ e, S) for e in E]
    P = np.zeros((r, r)); P[S['piR'], np.arange(r)] = 1                    # c -> coefficients of R(K)
    LC = [L @ P for L in LR]
    nI = len(E)
    print(f"ns={ns} N={N} orbits r={r} center dim={zc} irreps={nI} weights={w}  ({time.time()-t0:.1f}s)", flush=True)
    # SWAP blocks: sign / whether SWAP E_rho = +-E_rho
    pts = {}
    order = sorted(range(nI), key=lambda i: -w[i])
    subsetsC = []
    if nI <= maxsub:
        for m in range(nI+1):
            subsetsC += list(itertools.combinations(range(nI), m))
    else:
        for _ in range(NS):
            subsetsC.append(tuple(i for i in range(nI) if rng.random() < rng.random()))
    best = {}
    for ZC in subsetsC:
        if ZC:
            Mc = np.vstack([LC[i] for i in ZC])
            u_, s_, vt = np.linalg.svd(Mc)
            nb = vt[np.sum(s_ > 1e-9 * max(1, s_[0] if len(s_) else 1)):].T
        else:
            nb = np.eye(r)
        if nb.shape[1] == 0: continue
        # greedy Z_R
        for start in range(3):
            ordR = order if start == 0 else list(rng.permutation(nI))
            ZR = []; A_ = np.zeros((0, nb.shape[1])); b_ = np.zeros(0)
            for i in ordR:
                A2 = np.vstack([A_, LR[i] @ nb]); b2 = np.concatenate([b_, sR[i]])
                z, *_ = np.linalg.lstsq(A2, b2, rcond=None)
                if np.linalg.norm(A2 @ z - b2) < 1e-7 * (1 + np.linalg.norm(b2)):
                    ZR.append(i); A_, b_ = A2, b2
            z, *_ = np.linalg.lstsq(A_, b_, rcond=None) if len(b_) else (np.zeros(nb.shape[1]),)
            # generic point
            uu, ss, vv = np.linalg.svd(A_) if len(b_) else (None, np.zeros(0), np.zeros((0, nb.shape[1])))
            nk = vv[np.sum(ss > 1e-9):].T if len(b_) else np.eye(nb.shape[1])
            zg = z + (nk @ rng.standard_normal(nk.shape[1]) if nk.shape[1] else 0)
            c = nb @ zg
            K = mat(c, S)
            p = rk(SW - K); Rm = K.reshape(N, N, N, N).transpose(0, 2, 1, 3).reshape(D, D); q = rk(Rm)
            if q not in best or p < best[q][0]:
                best[q] = (p, ZC, tuple(ZR), c.tolist())
    # Pareto
    fr = []; bp = 10**9
    for q in sorted(best):
        if best[q][0] < bp:
            bp = best[q][0]; fr.append((best[q][0], q))
    print("frontier (p,q):", fr, f"({time.time()-t0:.1f}s)", flush=True)
    if out:
        json.dump({str(q): [int(best[q][0]), [int(i) for i in best[q][1]], [int(i) for i in best[q][2]], best[q][3]] for q in best}, open(out, 'w'))
    return fr

NS=3000
if __name__ == '__main__':
    if len(sys.argv)>2: NS=int(sys.argv[2])
    ns = list(map(int, sys.argv[1].split(',')))
    search(ns, out=f"wreath_{'_'.join(map(str,ns))}.json")
