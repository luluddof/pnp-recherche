"""Proof-size measurements for hard CNF families.
Measures (all are UPPER bounds on minimal proof size in the given system):
  * CDCL (Glucose3) conflicts and DRUP proof length (lines)  -> general resolution (RUP ~ resolution)
  * DPLL tree leaves (simple heuristic)                       -> tree-like resolution
  * Gaussian elimination row-ops for XOR families             -> 'resolution + XOR' / Extended Frege
"""
import itertools, random, sys, time, statistics
import networkx as nx
from pysat.solvers import Solver

sys.setrecursionlimit(100000)

# ---------------- generators ----------------
def php(n):  # n+1 pigeons, n holes
    v = lambda i, j: i * n + j + 1
    cls = [[v(i, j) for j in range(n)] for i in range(n + 1)]
    for j in range(n):
        for a, b in itertools.combinations(range(n + 1), 2):
            cls.append([-v(a, j), -v(b, j)])
    return cls

def xor_clauses(vars_, parity):  # sum vars = parity (mod 2)
    out = []
    k = len(vars_)
    for signs in itertools.product([0, 1], repeat=k):
        if sum(signs) % 2 != parity:  # forbid assignments with wrong parity
            # clause falsified exactly by assignment 'signs'
            out.append([(-v if s else v) for v, s in zip(vars_, signs)])
    return out

def tseitin(nv, d=4, seed=0):
    G = nx.random_regular_graph(d, nv, seed=seed)
    eid = {frozenset(e): i + 1 for i, e in enumerate(G.edges())}
    charge = {u: 0 for u in G}
    charge[0] = 1  # total charge odd -> unsat
    cls, rows = [], []
    for u in G:
        vs = [eid[frozenset((u, w))] for w in G[u]]
        cls += xor_clauses(vs, charge[u])
        rows.append((set(vs), charge[u]))
    return cls, rows, len(eid)

def rand3xor(n, m, seed):
    """Random 3-XOR system, m>n equations -> unsat w.h.p.; we force unsat by checking."""
    rng = random.Random(seed)
    rows = []
    for _ in range(m):
        vs = rng.sample(range(1, n + 1), 3)
        rows.append((set(vs), rng.randint(0, 1)))
    cls = []
    for vs, b in rows:
        cls += xor_clauses(sorted(vs), b)
    return cls, rows

def rand3cnf(n, ratio, seed):
    rng = random.Random(seed)
    return [[v * rng.choice([1, -1]) for v in rng.sample(range(1, n + 1), 3)]
            for _ in range(int(round(ratio * n)))]

# ---------------- measures ----------------
def cdcl(cls):
    with Solver(name='g3', bootstrap_with=cls, with_proof=True) as s:
        sat = s.solve()
        st = s.accum_stats()
        prf = s.get_proof() if not sat else None
    return sat, st['conflicts'], (len(prf) if prf is not None else None)

def dpll_leaves(cls, limit=2_000_000):
    """Count leaves of a DPLL tree (= leaves of a tree-like resolution refutation)."""
    cls = [frozenset(c) for c in cls]
    cnt = [0]

    def simplify(cl, lit):
        out = []
        for c in cl:
            if lit in c:
                continue
            if -lit in c:
                c = c - {-lit}
                if not c:
                    return None
            out.append(c)
        return out

    def rec(cl):
        if cnt[0] > limit:
            raise RuntimeError
        # unit propagation
        while True:
            unit = next((c for c in cl if len(c) == 1), None)
            if unit is None:
                break
            cl = simplify(cl, next(iter(unit)))
            if cl is None:
                cnt[0] += 1
                return False
        if not cl:
            return True
        # branch: most frequent var in shortest clauses
        m = min(len(c) for c in cl)
        freq = {}
        for c in cl:
            if len(c) == m:
                for l in c:
                    freq[abs(l)] = freq.get(abs(l), 0) + 1
        x = max(freq, key=freq.get)
        for lit in (x, -x):
            sub = simplify(cl, lit)
            if sub is None:
                cnt[0] += 1
                continue
            if rec(sub):
                return True
        return False

    try:
        sat = rec(cls)
    except RuntimeError:
        return None, '>%d' % limit
    return sat, cnt[0]

def gauss_ops(rows):
    """Gaussian elimination over GF(2); returns (#row additions to derive 0=1) or None if consistent."""
    rows = [(frozenset(s), b) for s, b in rows]
    ops = 0
    pivots = []  # (pivot var, row)
    for s, b in rows:
        s = set(s)
        for pv, (ps, pb) in pivots:
            if pv in s:
                s ^= ps
                b ^= pb
                ops += 1
        if not s:
            if b == 1:
                return ops
            continue
        pv = min(s)
        pivots.append((pv, (frozenset(s), b)))
    return None

# ---------------- experiments ----------------
def main():
    print("== PHP(n+1,n) ==")
    print("n vars clauses | CDCL_conflicts DRUP_lines | DPLL_tree_leaves")
    for n in range(4, 11):
        cls = php(n)
        t = time.time()
        sat, conf, plen = cdcl(cls)
        tl = dpll_leaves(cls, 300_000)[1] if n <= 8 else '-'
        print(n, n * (n + 1), len(cls), '|', conf, plen, '|', tl, ' (%.1fs)' % (time.time() - t), flush=True)

    print("\n== Tseitin, random 4-regular graph, odd charge ==")
    print("|V| edges clauses | CDCL_conflicts DRUP_lines | DPLL_leaves | GaussOps")
    for nv in [10, 14, 18, 22, 26, 30, 34, 38, 42]:
        res = []
        for seed in range(3):
            cls, rows, ne = tseitin(nv, 4, seed)
            sat, conf, plen = cdcl(cls)
            assert not sat
            tl = dpll_leaves(cls, 300_000)[1] if nv <= 22 else '-'
            res.append((conf, plen, tl, gauss_ops(rows)))
        med = lambda k: statistics.median([r[k] for r in res])
        print(nv, ne, len(cls), '|', med(0), med(1), '|', [r[2] for r in res], '|', med(3), flush=True)

    print("\n== random 3-XOR, m = 1.5n equations (unsat instances) ==")
    print("n | CDCL_conflicts DRUP_lines | GaussOps")
    for n in [20, 40, 60, 80, 100, 120, 140]:
        res = []
        seed = 0
        while len(res) < 3:
            cls, rows = rand3xor(n, int(1.5 * n), 1000 * n + seed)
            seed += 1
            g = gauss_ops(rows)
            if g is None:
                continue
            sat, conf, plen = cdcl(cls)
            assert not sat
            res.append((conf, plen, g))
        print(n, '|', statistics.median(r[0] for r in res), statistics.median(r[1] for r in res), '|',
              statistics.median(r[2] for r in res), flush=True)

    print("\n== random 3-CNF, ratio 4.26 (unsat instances only) and ratio 6 ==")
    print("n ratio | #unsat/tried | median CDCL_conflicts | median DRUP_lines")
    for ratio in [4.26, 6.0]:
        for n in [50, 100, 150, 200, 250]:
            res = []
            tried = 0
            while len(res) < 5 and tried < 30:
                cls = rand3cnf(n, ratio, 7 * n + tried)
                tried += 1
                sat, conf, plen = cdcl(cls)
                if not sat:
                    res.append((conf, plen))
            if res:
                print(n, ratio, '|', '%d/%d' % (len(res), tried), '|',
                      statistics.median(r[0] for r in res), statistics.median(r[1] for r in res), flush=True)

if __name__ == '__main__':
    main()
