"""BLP, AIP (over Z, exact), AIP over GF(2), and BLP+AIP for Boolean CSP instances.
Ground truth by brute force. Compares on templates 2-SAT, XORSAT, Horn-3SAT, 3-SAT, 1-in-3 SAT.
"""
import itertools, random, sys
import numpy as np
from scipy.optimize import linprog

# ---------------- templates: list of relations (as sets of tuples) -------------
def rel(arity, pred):
    return frozenset(t for t in itertools.product((0, 1), repeat=arity) if pred(t))

def clause(signs):  # signs: tuple of 0/1, literal i true iff t[i]==signs[i]
    return rel(len(signs), lambda t: any(t[i] == signs[i] for i in range(len(signs))))

TEMPLATES = {
    '2-SAT': [clause(s) for s in itertools.product((0, 1), repeat=2)],
    'XORSAT': [rel(3, lambda t: sum(t) % 2 == 0), rel(3, lambda t: sum(t) % 2 == 1),
               rel(2, lambda t: t[0] != t[1])],
    'Horn': [clause((1,)), clause((0,)), clause((0, 0, 1)), clause((0, 0, 0)), clause((0, 1))],
    '3-SAT': [clause(s) for s in itertools.product((0, 1), repeat=3)],
    '1-in-3': [rel(3, lambda t: sum(t) == 1)],
    # pure-negative/positive 3-SAT style "NAE" for extra contrast
    'NAE-3': [rel(3, lambda t: 0 < sum(t) < 3)],
}

# ---------------- instance = (n, [(scope, relation)]) --------------------------
def brute(n, cons):
    for a in itertools.product((0, 1), repeat=n):
        if all(tuple(a[v] for v in sc) in R for sc, R in cons):
            return True
    return False

def build_system(n, cons):
    """Variables: p[x,a] for x<n, a in {0,1}; q[c,t] for t in R_c.
    Equations: sum_a p[x,a]=1; sum_t q[c,t]=1; sum_{t:t_i=a} q[c,t] = p[s_i,a]."""
    idx = {}
    for x in range(n):
        for a in (0, 1):
            idx[('p', x, a)] = len(idx)
    for c, (sc, R) in enumerate(cons):
        for t in sorted(R):
            idx[('q', c, t)] = len(idx)
    rows, rhs = [], []
    N = len(idx)
    for x in range(n):
        r = [0] * N; r[idx[('p', x, 0)]] = 1; r[idx[('p', x, 1)]] = 1
        rows.append(r); rhs.append(1)
    for c, (sc, R) in enumerate(cons):
        r = [0] * N
        for t in R: r[idx[('q', c, t)]] = 1
        rows.append(r); rhs.append(1)
        for i, v in enumerate(sc):
            for a in (0, 1):
                r = [0] * N
                for t in R:
                    if t[i] == a: r[idx[('q', c, t)]] += 1
                r[idx[('p', v, a)]] -= 1
                rows.append(r); rhs.append(0)
    return idx, rows, rhs

def blp(n, cons):
    idx, A, b = build_system(n, cons)
    res = linprog(np.zeros(len(idx)), A_eq=np.array(A, float), b_eq=np.array(b, float),
                  bounds=[(0, None)] * len(idx), method='highs')
    return res.status == 0

def blp_support(n, cons):
    """Set of variables that are positive in some BLP solution (= support of relative-interior pt)."""
    idx, A, b = build_system(n, cons)
    A = np.array(A, float); b = np.array(b, float)
    N = len(idx); supp = set()
    for j in range(N):
        if j in supp: continue
        c = np.zeros(N); c[j] = -1
        res = linprog(c, A_eq=A, b_eq=b, bounds=[(0, 1)] * N, method='highs')
        if res.status != 0: return None
        if -res.fun > 1e-9:
            x = res.x
            supp |= {k for k in range(N) if x[k] > 1e-9}
    return supp

# exact integer solvability of A x = b via column-HNF (unimodular column ops)
def int_solvable(A, b):
    m = len(A)
    if m == 0: return True
    N = len(A[0])
    cols = [[A[i][j] for i in range(m)] for j in range(N)]
    b = list(b)
    piv = 0
    for i in range(m):
        # gcd-reduce entries in row i among columns piv..end
        while True:
            nz = [j for j in range(piv, len(cols)) if cols[j][i] != 0]
            if not nz: break
            jmin = min(nz, key=lambda j: abs(cols[j][i]))
            cols[piv], cols[jmin] = cols[jmin], cols[piv]
            done = True
            for j in range(piv + 1, len(cols)):
                if cols[j][i] != 0:
                    qq = cols[j][i] // cols[piv][i]
                    cols[j] = [cols[j][k] - qq * cols[piv][k] for k in range(m)]
                    if cols[j][i] != 0: done = False
            if done: break
        if any(cols[j][i] != 0 for j in range(piv, len(cols))):
            # pivot column piv; eliminate b's row i
            d = cols[piv][i]
            if b[i] % d != 0: return False
            qq = b[i] // d
            b = [b[k] - qq * cols[piv][k] for k in range(m)]
            piv += 1
        else:
            if b[i] != 0: return False
    return all(v == 0 for v in b)

def aip(n, cons, restrict=None):
    idx, A, b = build_system(n, cons)
    if restrict is not None:
        keep = sorted(restrict)
        A = [[r[j] for j in keep] for r in A]
    return int_solvable(A, b)

def aip_gf2(n, cons):
    idx, A, b = build_system(n, cons)
    M = [[x % 2 for x in r] + [y % 2] for r, y in zip(A, b)]
    N = len(idx); r = 0
    for c in range(N):
        p = next((i for i in range(r, len(M)) if M[i][c]), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        for i in range(len(M)):
            if i != r and M[i][c]:
                M[i] = [x ^ y for x, y in zip(M[i], M[r])]
        r += 1
    return not any(all(x == 0 for x in row[:-1]) and row[-1] for row in M)

def blp_aip(n, cons):
    s = blp_support(n, cons)
    if s is None: return False
    return aip(n, cons, restrict=s)

# ---------------- instance generation -------------------------------------------
def rand_instance(T, n, m, rng):
    cons = []
    for _ in range(m):
        R = rng.choice(T)
        k = len(next(iter(R)))
        sc = tuple(rng.sample(range(n), k))
        cons.append((sc, R))
    return cons

def run(trials=400, seed=1):
    rng = random.Random(seed)
    algs = {'BLP': blp, 'AIP_Z': aip, 'AIP_GF2': aip_gf2, 'BLP+AIP': blp_aip}
    out = {}
    for name, T in TEMPLATES.items():
        stats = {a: [0, 0] for a in algs}  # [correct, false-accept]
        nsat = 0
        for _ in range(trials):
            n = rng.randint(3, 7)
            m = rng.randint(n // 2 + 1, 3 * n)
            cons = rand_instance(T, n, m, rng)
            truth = brute(n, cons); nsat += truth
            for a, f in algs.items():
                r = f(n, cons)
                if r == truth: stats[a][0] += 1
                elif r and not truth: stats[a][1] += 1
                else: raise RuntimeError(f'relaxation rejected satisfiable instance?! {a} {name}')
        out[name] = (stats, nsat)
        print(name, f'sat={nsat}/{trials}',
              ' '.join(f'{a}:{c}/{trials}(FA={fa})' for a, (c, fa) in stats.items()), flush=True)
    return out

if __name__ == '__main__':
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 300)
