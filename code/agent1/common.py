"""Instance generators, ground truth (CaDiCaL via python-sat), unit propagation."""
import random, itertools
from pysat.solvers import Cadical153


def rand3sat(n, m, rng):
    cls = []
    for _ in range(m):
        vs = rng.sample(range(1, n + 1), 3)
        cls.append([v if rng.random() < .5 else -v for v in vs])
    return n, cls


def xor_clauses(vs, parity):
    """CNF for x_vs[0] xor ... = parity (all 2^(k-1) clauses)."""
    out = []
    k = len(vs)
    for signs in itertools.product([0, 1], repeat=k):
        # forbid assignment a = signs if its parity != parity
        if sum(signs) % 2 != parity:
            out.append([(-v if s else v) for v, s in zip(vs, signs)])
    return out


def rand3xor(n, m, rng):
    cls = []
    for _ in range(m):
        vs = rng.sample(range(1, n + 1), 3)
        cls += xor_clauses(vs, rng.randrange(2))
    return n, cls


def random_cubic_graph(nv, rng):
    while True:
        pts = [v for v in range(nv) for _ in range(3)]
        rng.shuffle(pts)
        edges = [(pts[2 * i], pts[2 * i + 1]) for i in range(len(pts) // 2)]
        if all(a != b for a, b in edges) and len(set(tuple(sorted(e)) for e in edges)) == len(edges):
            return edges


def tseitin(nv, rng, odd=True):
    """Tseitin on random 3-regular multigraph-free graph; total charge odd => UNSAT."""
    edges = random_cubic_graph(nv, rng)
    inc = [[] for _ in range(nv)]
    for i, (a, b) in enumerate(edges):
        inc[a].append(i + 1); inc[b].append(i + 1)
    charge = [rng.randrange(2) for _ in range(nv)]
    if (sum(charge) % 2 == 1) != odd:
        charge[0] ^= 1
    cls = []
    for v in range(nv):
        cls += xor_clauses(inc[v], charge[v])
    return len(edges), cls


def php(h):
    """h+1 pigeons into h holes. var p(i,j) = i*h+j+1."""
    P = h + 1
    var = lambda i, j: i * h + j + 1
    cls = [[var(i, j) for j in range(h)] for i in range(P)]
    for j in range(h):
        for i in range(P):
            for k in range(i + 1, P):
                cls.append([-var(i, j), -var(k, j)])
    return P * h, cls


def truth(n, cls):
    with Cadical153(bootstrap_with=cls) as s:
        return s.solve()


def check(cls, assign):
    """assign: dict var->bool or list index var."""
    return all(any((assign[abs(l)] if l > 0 else not assign[abs(l)]) for l in c) for c in cls)


def unit_prop(cls, assign):
    """assign: dict var->bool (modified copy returned) or None on conflict."""
    a = dict(assign)
    changed = True
    while changed:
        changed = False
        for c in cls:
            unk = None; nunk = 0; sat = False
            for l in c:
                v = abs(l)
                if v in a:
                    if a[v] == (l > 0):
                        sat = True; break
                else:
                    nunk += 1; unk = l
            if sat:
                continue
            if nunk == 0:
                return None
            if nunk == 1:
                a[abs(unk)] = unk > 0; changed = True
    return a
