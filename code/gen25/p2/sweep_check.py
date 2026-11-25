# Check of the sweep lemma (level 0 of Theorem A by cut-width sweep) against brute force.
# (a) random integer tables on a (K,lam)-clique: sweep value == brute-force value, and op counts.
# (b) end-to-end on random small 3-CNFs: clauses -> K groups (Lemma 1 / Lemma 8 construction of the
#     (K,L)-clique tables), sweep over the clique == brute-force #SAT over the occurring variables.
import itertools, random, sys
from collections import defaultdict

def edges_of(K): return [(a, b) for a in range(K) for b in range(a+1, K)]

def brute_clique(K, lam, W):
    """W[w](dict: neighbour -> string) ; sum over one string per edge of prod_w W_w."""
    E = edges_of(K); tot = 0; ops = 0
    for strs in itertools.product(range(2**lam), repeat=len(E)):
        s = dict(zip(E, strs)); prod = 1
        for w in range(K):
            prod *= W(w, {v: s[(min(w, v), max(w, v))] for v in range(K) if v != w})
            if prod == 0: break
        tot += prod; ops += 1
    return tot, ops

def sweep_clique(K, lam, W):
    """Sweep nodes 0..K-1. F_j is a table on the strings of the cut edges {0..j-1} x {j..K-1}.
       Step j enumerates the strings on cut_j (= edges (w,v), w<j<=v) plus the new edges (j,v), v>j."""
    F = {(): 1}; cut = []   # cut: ordered list of edges (w,v), w<j<=v
    ops = 0; maxcut = 0
    for j in range(K):
        new = [(j, v) for v in range(j+1, K)]
        inc = [e for e in cut if e[1] == j]          # edges into j: summed out now
        keep = [e for e in cut if e[1] != j]
        ncut = keep + new
        maxcut = max(maxcut, len(cut) + len(new))
        G = defaultdict(int)
        for key, val in F.items():
            s = dict(zip(cut, key))
            for nv in itertools.product(range(2**lam), repeat=len(new)):
                ops += 1
                sn = dict(zip(new, nv))
                arg = {}
                for (w, v) in inc: arg[w] = s[(w, v)]
                for (w, v) in new: arg[v] = sn[(w, v)]
                x = W(j, arg)
                if x:
                    G[tuple(s[e] for e in keep) + nv] += val * x
        F = dict(G); cut = ncut
    return F.get((), 0), ops, maxcut

def part_a(trials=30, seed=1):
    rng = random.Random(seed); ok = True
    for t in range(trials):
        K = rng.choice([3, 4, 5]); lam = 1 if K == 5 else rng.choice([1, 2])
        tables = [dict() for _ in range(K)]
        def W(w, arg, tables=tables):
            key = tuple(arg[v] for v in sorted(arg))
            if key not in tables[w]: tables[w][key] = rng.randint(-3, 3)
            return tables[w][key]
        # fix tables first (deterministic contents independent of query order)
        for w in range(K):
            for key in itertools.product(range(2**lam), repeat=K-1): tables[w][key] = rng.randint(-3, 3)
        zb, ob = brute_clique(K, lam, W); zs, os_, mc = sweep_clique(K, lam, W)
        pred = ((K+1)**2)//4 - 1
        ok &= (zb == zs) and (mc == pred)
        if t < 8 or zb != zs:
            print(f"  (a) K={K} lam={lam}: brute={zb} sweep={zs} equal={zb==zs}  brute terms={ob}  sweep terms={os_}"
                  f"  max edges enumerated in a step={mc} (formula floor((K+1)^2/4)-1={pred}; brute K(K-1)/2={K*(K-1)//2})")
    print(f" (a) all {trials} random-table instances agree: {ok}")
    return ok

def rand3cnf(n, m, rng):
    F = []
    while len(F) < m:
        vs = rng.sample(range(n), 3)
        c = frozenset((v, rng.random() < 0.5) for v in vs)
        if c not in F: F.append(c)
    return [sorted(c) for c in F]

def count_brute(F, n):
    occ = sorted({v for c in F for v, _ in c}); cnt = 0
    for bits in itertools.product([0, 1], repeat=len(occ)):
        a = dict(zip(occ, bits))
        if all(any(a[v] == (1 if pos else 0) for v, pos in c) for c in F): cnt += 1
    return cnt

def cnf_clique(F, K, rng):
    """Network of F (Section 6 of the paper): one bit per occurrence, chain tests per variable;
       clauses -> K groups (random). Returns K, L and the tables of the (K,L)-clique (Lemma 1 / eq. (26))."""
    m = len(F); group = [i % K for i in range(m)]; rng.shuffle(group)
    occs = [(ci, v, pos) for ci, c in enumerate(F) for v, pos in c]      # occurrence o = index
    byvar = defaultdict(list)
    for o, (ci, v, pos) in enumerate(occs): byvar[v].append(o)
    tests = [(lst[i], lst[i+1]) for lst in byvar.values() for i in range(len(lst)-1)]
    g = lambda o: group[occs[o][0]]
    cross = defaultdict(list)   # (gamma,gamma') gamma<gamma' -> ordered list of tests (block positions)
    internal = defaultdict(list)
    for (o1, o2) in tests:
        a, b = g(o1), g(o2)
        if a == b: internal[a].append((o1, o2))
        else:
            key = (min(a, b), max(a, b)); cross[key].append((o1, o2))
    L = max([len(v) for v in cross.values()] + [1])
    myocc = defaultdict(list)
    for o in range(len(occs)): myocc[g(o)].append(o)
    myclauses = defaultdict(list)
    for ci in range(m): myclauses[group[ci]].append(ci)
    cache = {}
    def W(gam, arg):
        # arg[gamma'] = integer whose bits are the L positions of the edge string {gam,gamma'}
        key = (gam, tuple(sorted(arg.items())))
        if key in cache: return cache[key]
        forced = {}
        for gp, s in arg.items():
            lst = cross.get((min(gam, gp), max(gam, gp)), [])
            if s >> len(lst): cache[key] = 0; return 0        # appended bits must be 0
            for pos, (o1, o2) in enumerate(lst):
                o = o1 if g(o1) == gam else o2
                b = (s >> pos) & 1
                if forced.get(o, b) != b: cache[key] = 0; return 0   # (cannot happen: one test per occurrence pair)
                forced[o] = b
        free = [o for o in myocc[gam] if o not in forced]; tot = 0
        for bits in itertools.product([0, 1], repeat=len(free)):
            z = dict(forced); z.update(zip(free, bits))
            if any(z[o1] != z[o2] for (o1, o2) in internal[gam]): continue
            if all(any(z[o] == (1 if occs[o][2] else 0) for o in myocc[gam] if occs[o][0] == ci) for ci in myclauses[gam]):
                tot += 1
        cache[key] = tot; return tot
    return L, W

def part_b(trials=40, seed=7):
    rng = random.Random(seed); ok = True; shown = 0
    for t in range(trials):
        n = rng.randint(5, 9); m = rng.randint(4, 9); K = rng.choice([3, 4])
        F = rand3cnf(n, m, rng)
        L, W = cnf_clique(F, K, rng)
        if 2**(L*((K+1)**2//4 - 1)) > 2**16: continue
        zs, os_, mc = sweep_clique(K, L, W); zb = count_brute(F, n)
        ok &= (zs == zb); shown += 1
        if shown <= 8 or zs != zb:
            print(f"  (b) n={n} m={m} K={K} L={L}: #SAT brute={zb} sweep={zs} equal={zs==zb}  sweep terms={os_}"
                  f"  vs brute-force clique terms 2^(L K(K-1)/2)={2**(L*K*(K-1)//2)}")
    print(f" (b) all {shown} random 3-CNF instances agree: {ok}")
    return ok

if __name__ == "__main__":
    a = part_a(); b = part_b()
    print("ALL OK" if a and b else "MISMATCH")
