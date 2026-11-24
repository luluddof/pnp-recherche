#!/usr/bin/env python3
"""G18i experiment: how small is the unsatisfiable part of a random 3-CNF at density D?
A size-oblivious refuter (the 2026 counting algorithms: cost depends only on the number of
clauses m', resp. n'+m') only needs ANY unsatisfiable sub-formula.  We measure
  (a) prefix core : smallest p such that the clauses inside variables 1..p are UNSAT (no search)
  (b) peel core   : same along a min-degree peeling order (greedy densest sub-hypergraph)
  (c) MUS         : deletion-based clause-minimal core inside (b) (upper bound on the min core)
and the CDCL (resolution) effort on full formula / cores, as a contrast.
"""
import sys, random, time, json, math
from pysat.solvers import Solver

def gen(n, D, seed):
    rng = random.Random(seed)
    m = int(round(D * n)); cls = []
    for _ in range(m):
        vs = rng.sample(range(1, n + 1), 3)
        cls.append([v if rng.random() < 0.5 else -v for v in vs])
    return cls

def solve(cls, name='cadical153'):
    with Solver(name=name, bootstrap_with=cls) as s:
        r = s.solve()
        st = s.accum_stats() or {}
    return r, st.get('conflicts', 0)

def smallest_unsat_prefix(cls, order):
    """order: list of variables; nested sets S_p = first p of order. Returns (p, induced clauses)."""
    pos = {v: i for i, v in enumerate(order)}
    key = [max(pos[abs(l)] for l in c) for c in cls]   # clause is inside S_p iff key < p
    def induced(p): return [c for c, k in zip(cls, key) if k < p]
    lo, hi = 0, len(order)                      # induced(lo) SAT, induced(hi) UNSAT
    if solve(induced(hi))[0]: return None, None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if solve(induced(mid))[0]: lo = mid
        else: hi = mid
    return hi, induced(hi)

def peel_order(n, cls):
    """min-degree peeling; returns variables in REVERSE removal order (last removed first)."""
    occ = {v: set() for v in range(1, n + 1)}
    for i, c in enumerate(cls):
        for l in c: occ[abs(l)].add(i)
    deg = {v: len(occ[v]) for v in occ}; alive = [True] * len(cls)
    import heapq
    h = [(d, v) for v, d in deg.items()]; heapq.heapify(h); removed = set(); out = []
    while h:
        d, v = heapq.heappop(h)
        if v in removed or d != deg[v]: continue
        removed.add(v); out.append(v)
        for i in occ[v]:
            if alive[i]:
                alive[i] = False
                for l in cls[i]:
                    u = abs(l)
                    if u not in removed:
                        deg[u] -= 1; heapq.heappush(h, (deg[u], u))
    return out[::-1]

def mus(n, cls, budget=40.0):
    """deletion-based MUS with core refinement; the working set is UNSAT at all times."""
    t0 = time.time(); m = len(cls)
    s = Solver(name='glucose4')
    for i, c in enumerate(cls): s.add_clause(c + [-(n + 1 + i)])
    work = list(range(m)); crit = set(); done = True
    assert not s.solve(assumptions=[n + 1 + i for i in work])
    core = s.get_core(); work = [i for i in work if (n + 1 + i) in set(core)]
    idx = 0
    while idx < len(work):
        if time.time() - t0 > budget: done = False; break
        i = work[idx]
        if i in crit: idx += 1; continue
        trial = [j for j in work if j != i]
        if s.solve(assumptions=[n + 1 + j for j in trial]):
            crit.add(i); idx += 1
        else:
            core = set(s.get_core()); work = [j for j in trial if (n + 1 + j) in core]; idx = 0
            work.sort(key=lambda j: j in crit)       # untested first
    s.delete()
    sub = [cls[i] for i in work]
    return sub, done

def nvars(cls): return len({abs(l) for c in cls for l in c})

def one(n, D, seed, do_mus=True):
    cls = gen(n, D, seed); r = {'n': n, 'D': D, 'seed': seed, 'm': len(cls)}
    sat, conf = solve(cls); r['full_unsat'] = (not sat); r['full_conf'] = conf
    if sat: return r
    p, pc = smallest_unsat_prefix(cls, list(range(1, n + 1)))
    r.update(pre_n=p, pre_m=len(pc), pre_conf=solve(pc)[1])
    q, qc = smallest_unsat_prefix(cls, peel_order(n, cls))
    r.update(peel_n=q, peel_m=len(qc), peel_conf=solve(qc)[1])
    if do_mus:
        mc, done = mus(n, qc)
        r.update(mus_m=len(mc), mus_n=nvars(mc), mus_done=done, mus_conf=solve(mc)[1])
    return r

if __name__ == '__main__':
    n = int(sys.argv[1]); Ds = [float(x) for x in sys.argv[2].split(',')]
    seeds = range(int(sys.argv[3]), int(sys.argv[4])); out = sys.argv[5]
    with open(out, 'a') as f:
        for D in Ds:
            for sd in seeds:
                t = time.time(); r = one(n, D, 1000 * sd + int(D * 10)); r['t'] = round(time.time() - t, 2)
                f.write(json.dumps(r) + '\n'); f.flush()
