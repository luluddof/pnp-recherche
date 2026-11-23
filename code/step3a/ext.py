#!/usr/bin/env python3
"""Step 3a: do extension variables (poly-time preprocessing) change the growth
rate of CaDiCaL refutations of UNSAT random 3-CNF?

Variants (all add definitions = Extended Resolution extension rule):
  base    : plain formula
  bva     : SimpleBVA (Manthey et al. 2012), only reductions >0
  bva0    : relaxed BVA also accepting 0-reduction replacements (more defs)
  xor     : y = a XOR b for the n most co-occurring variable pairs (defs only)
  orsub   : y <-> (l1 v l2) for most co-occurring literal pairs, substituted
            into every clause containing both literals (shortens clauses)
  and     : y <-> (a AND b) for most co-occurring var pairs (defs only)
  rand    : control: n random pairs, random op XOR/AND, defs only
Measure: CaDiCaL 1.9.5 conflicts and DRUP proof lines (+ #definition clauses).
"""
import sys, random, itertools, json, time, os
from collections import Counter, defaultdict
from multiprocessing import Pool
from pysat.solvers import Solver

DENS = 4.26


def gen(n, seed, d=DENS):
    rng = random.Random(seed * 1000003 + n)
    m = round(d * n)
    cls = []
    for _ in range(m):
        vs = rng.sample(range(1, n + 1), 3)
        cls.append([v if rng.random() < 0.5 else -v for v in vs])
    return cls


def solve(n, cls, proof=False):
    with Solver(name='cadical195', bootstrap_with=cls, with_proof=proof) as s:
        r = s.solve()
        st = s.accum_stats()
        pl = len(s.get_proof()) if (proof and not r) else None
    return r, st['conflicts'], pl


# ---------------- extension variants ----------------
def defs_xor(y, a, b):
    return [[-y, a, b], [-y, -a, -b], [y, -a, b], [y, a, -b]]


def defs_and(y, a, b):  # y <-> a & b  (a,b literals)
    return [[-y, a], [-y, b], [y, -a, -b]]


def defs_or(y, a, b):  # y <-> a v b
    return [[y, -a], [y, -b], [-y, a, b]]


def var_pairs(cls):
    c = Counter()
    for C in cls:
        for a, b in itertools.combinations(sorted(abs(l) for l in C), 2):
            c[(a, b)] += 1
    return c


def ext_xor(n, cls, k):
    pairs = [p for p, _ in var_pairs(cls).most_common(k)]
    out, y = list(cls), n
    for a, b in pairs:
        y += 1
        out += defs_xor(y, a, b)
    return y, out, 4 * len(pairs)


def ext_and(n, cls, k):
    pairs = [p for p, _ in var_pairs(cls).most_common(k)]
    out, y = list(cls), n
    for a, b in pairs:
        y += 1
        out += defs_and(y, a, b)
    return y, out, 3 * len(pairs)


def ext_rand(n, cls, k, seed):
    rng = random.Random(seed + 77)
    out, y, nd = list(cls), n, 0
    for _ in range(k):
        a, b = rng.sample(range(1, n + 1), 2)
        a = a if rng.random() < .5 else -a
        b = b if rng.random() < .5 else -b
        y += 1
        d = defs_xor(y, a, b) if rng.random() < .5 else defs_and(y, a, b)
        out += d
        nd += len(d)
    return y, out, nd


def ext_orsub(n, cls, k):
    """Greedy: repeatedly take most frequent literal pair, define y<->(a v b),
    substitute into all clauses containing both a and b."""
    cls = [list(C) for C in cls]
    y, defs = n, []
    for _ in range(k):
        c = Counter()
        for C in cls:
            for a, b in itertools.combinations(sorted(C), 2):
                c[(a, b)] += 1
        if not c:
            break
        (a, b), f = c.most_common(1)[0]
        if f < 2:
            break
        y += 1
        defs += defs_or(y, a, b)
        for C in cls:
            if a in C and b in C:
                C.remove(a); C.remove(b); C.append(y)
    return y, cls + defs, len(defs)


def ext_bva(n, cls, allow_zero=False, maxdefs=None):
    """SimpleBVA (Manthey, Heule, Biere 2012)."""
    F = set(frozenset(C) for C in cls)
    occ = defaultdict(set)
    for C in F:
        for l in C:
            occ[l].add(C)
    y = n
    nd = 0
    queue = sorted(occ, key=lambda l: -len(occ[l]))
    tried = set()
    while queue:
        if maxdefs is not None and nd >= maxdefs:
            break
        l = queue.pop(0)
        if not occ[l]:
            continue
        Mlit = [l]
        Mcls = set(occ[l])
        while True:
            P = []  # (lit, clauseC in Mcls)
            for C in Mcls:
                rest = C - {l}
                lmin = min(rest, key=lambda z: len(occ[z])) if rest else None
                if lmin is None:
                    continue
                for D in occ[lmin]:
                    if len(D) == len(C):
                        diff = D - C
                        if len(diff) == 1 and (C - D) == {l}:
                            P.append((next(iter(diff)), C))
            cnt = Counter(lp for lp, _ in P if lp not in Mlit)
            if not cnt:
                break
            lmax, _ = cnt.most_common(1)[0]
            newMcls = set(C for lp, C in P if lp == lmax)
            old = len(Mlit) * len(Mcls) - len(Mlit) - len(Mcls)
            new = (len(Mlit) + 1) * len(newMcls) - (len(Mlit) + 1) - len(newMcls)
            if new > old:
                Mlit.append(lmax); Mcls = newMcls
            else:
                break
        red = len(Mlit) * len(Mcls) - len(Mlit) - len(Mcls)
        if len(Mlit) == 1 or (red < 0) or (red == 0 and not allow_zero):
            continue
        key = (tuple(sorted(Mlit)), frozenset(Mcls))
        if key in tried:
            continue
        tried.add(key)
        y += 1
        nd += 1
        for li in Mlit:
            occ  # replace
            for C in Mcls:
                D = (C - {l}) | {li}
                if D in F:
                    F.discard(D)
                    for z in D:
                        occ[z].discard(D)
            N = frozenset([li, y])
            F.add(N)
            for z in N:
                occ[z].add(N)
        for C in Mcls:
            N = (C - {l}) | {-y}
            F.add(N)
            for z in N:
                occ[z].add(N)
        queue.insert(0, l)
        for li in Mlit:
            if li not in queue:
                queue.append(li)
    return y, [list(C) for C in F], nd


def variants(n, cls, seed):
    yield 'base', (n, cls, 0)
    yield 'bva', ext_bva(n, cls)
    yield 'bva0', ext_bva(n, cls, allow_zero=True)
    yield 'xor', ext_xor(n, cls, n)
    yield 'and', ext_and(n, cls, n)
    yield 'orsub', ext_orsub(n, cls, n)
    yield 'rand', ext_rand(n, cls, n, seed)


def job(args):
    n, seed = args
    cls = gen(n, seed)
    r, c, pl = solve(n, cls, proof=True)
    if r:
        return None
    res = {'n': n, 'seed': seed}
    for name, (nv, ecls, nd) in variants(n, cls, seed):
        t = time.time()
        r2, c2, pl2 = solve(nv, ecls, proof=True)
        assert not r2, name
        res[name] = {'conf': c2, 'proof': pl2, 'defs': nd, 'nv': nv,
                     't': round(time.time() - t, 2)}
    return res


if __name__ == '__main__':
    ns = [int(x) for x in sys.argv[1].split(',')]
    per = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    out = sys.argv[3] if len(sys.argv) > 3 else 'results.jsonl'
    with Pool(4) as pool, open(out, 'a') as f:
        for n in ns:
            got, seed = 0, 0
            while got < per:
                batch = [(n, s) for s in range(seed, seed + 2 * per)]
                seed += 2 * per
                for res in pool.imap(job, batch):
                    if res and got < per:
                        got += 1
                        f.write(json.dumps(res) + '\n'); f.flush()
            print('done n=', n, flush=True)
