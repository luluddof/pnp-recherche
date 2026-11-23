import random, sys
sys.setrecursionlimit(10000)

def gen(n, m, rng):
    return [tuple(rng.choice([1,-1])*v for v in rng.sample(range(1,n+1),3)) for _ in range(m)]

def propagate(cls, assign):
    # unit propagation; returns new assign or None on conflict
    a = dict(assign)
    changed = True
    while changed:
        changed = False
        for c in cls:
            unassigned, sat = [], False
            for l in c:
                v = a.get(abs(l))
                if v is None: unassigned.append(l)
                elif v == (l > 0): sat = True; break
            if sat: continue
            if not unassigned: return None
            if len(unassigned) == 1:
                l = unassigned[0]; a[abs(l)] = l > 0; changed = True
    return a

# Candidate polynomial algorithm: failed-literal probing + greedy choice, never backtracks.
def candidate(cls, n):
    a = propagate(cls, {})
    if a is None: return "UNSAT"
    while len(a) < n:
        progress = False
        for v in range(1, n+1):
            if v in a: continue
            t = propagate(cls, {**a, v: True}); f = propagate(cls, {**a, v: False})
            if t is None and f is None: return "UNSAT"
            if t is None: a = f; progress = True
            elif f is None: a = t; progress = True
        if not progress:
            v = next(v for v in range(1, n+1) if v not in a)
            # greedy: pick value satisfying most clauses
            score = sum((1 if l == v else -1 if l == -v else 0) for c in cls for l in c)
            a = propagate(cls, {**a, v: score >= 0})
            if a is None: return "UNKNOWN"   # greedy failed, no backtracking allowed
    return "SAT"

nodes = 0
def dpll(cls, n, a):
    global nodes; nodes += 1
    a = propagate(cls, a)
    if a is None: return False
    if len(a) == n: return True
    v = next(v for v in range(1, n+1) if v not in a)
    return dpll(cls, n, {**a, v: True}) or dpll(cls, n, {**a, v: False})

rng = random.Random(1)
print(" n | instances | SAT reels | candidat correct | noeuds DPLL moyens")
for n in [20, 30, 40, 50, 60]:
    m = round(4.26*n); ok = sat = 0; tot = 0; T = 30
    for _ in range(T):
        cls = gen(n, m, rng)
        nodes = 0; truth = dpll(cls, n, {}); tot += nodes
        r = candidate(cls, n); sat += truth
        ok += (r == "SAT" and truth) or (r == "UNSAT" and not truth)
    print(f"{n:2} | {T:9} | {sat:9} | {ok:16} | {tot/T:10.0f}")
