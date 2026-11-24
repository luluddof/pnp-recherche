# Theorem C (G7b): tau <= W^a + s - 1 + |Z|, and the "moreover" clause (tau + 2 - s - |Z|).
import random, itertools, sys
sys.path.insert(0, '.')
from c3_H_rado import rank_mod, formula, P
from functools import lru_cache
def wtw(nv, adj, w):
    """min over elimination orders of max weight of {v} + higher neighbourhood  (= min over TDs of max bag weight)."""
    if nv == 0: return 0
    nb = [sum(1 << u for u in adj[v]) for v in range(nv)]
    def Q(Smask, v):  # vertices outside S+{v} reachable from v through S
        seen = 1 << v; stack = [v]; out = 0
        while stack:
            x = stack.pop(); m = nb[x] & ~seen
            while m:
                y = (m & -m).bit_length() - 1; m &= m - 1; seen |= 1 << y
                if Smask >> y & 1: stack.append(y)
                else: out |= 1 << y
        return out
    @lru_cache(None)
    def T(Smask):
        if Smask == 0: return 0
        best = 10**9; m = Smask
        while m:
            v = (m & -m).bit_length() - 1; m &= m - 1
            q = Q(Smask & ~(1 << v), v)
            bw = w[v] + sum(w[u] for u in range(nv) if q >> u & 1)
            if bw >= best: continue
            best = min(best, max(bw, T(Smask & ~(1 << v))))
        return best
    return T((1 << nv) - 1)
def tw(nv, adj): return wtw(nv, adj, [1]*nv) - 1
# ---- counterexample to the "moreover" clause -------------------------------------------------
s = 3; n = 2*s
cl = list(itertools.combinations(range(n), 3))
lab = [0]*s + [1]*s            # group 0 unmeasured (R=0), group 1 fully measured; alpha(e)=0 for all e
# P_{e,v}(t) = e_v : restriction to group 1 is e_v or 0, in R_1 = Q^{S_1}; own coefficients all 1 => Z = {}
padj = {v: set(range(n)) - {v} for v in range(n)}
tau = tw(n, padj)
# I_t: edge 0-1; measured = {1}; torso = single vertex {1}; its only bag has sum a = s
Wa = s
print(f"counterexample n={n}, s={s}: tau={tau}, W^a={Wa}, |Z|=0")
print(f"   main bound  tau <= W^a + s - 1 + |Z| : {tau} <= {Wa+s-1}  ->", tau <= Wa+s-1)
print(f"   'moreover'  max bag sum a >= tau + 2 - s - |Z| : {Wa} >= {tau+2-s}  ->", Wa >= tau+2-s, "  (corrected form tau+1-s-|Z|:", Wa >= tau+1-s, ")")
# ---- random per-point structures ---------------------------------------------------------------
rng = random.Random(3); cnt = 0; minslack = 99; mor_fail = 0
for it in range(150):
    n = rng.choice([9, 10, 12]); g = rng.choice([3, 4, 5, 6])
    cl = formula(n, rng)
    lab = [rng.randrange(g) for _ in range(n)]
    S = [[v for v in range(n) if lab[v] == h] for h in range(g)]
    a = [len(x) for x in S]
    meas = [a[h] > 0 and rng.random() < .6 for h in range(g)]
    alpha = []
    for e in cl:
        um = [lab[v] for v in e if not meas[lab[v]]]
        alpha.append(rng.choice(um) if um and rng.random() < .9 else lab[rng.choice(e)])
    adj = {h: set() for h in range(g)}; Z = set(); F2 = []
    for e, vs in enumerate(cl):
        clique = {alpha[e]}; keep = []
        for v in vs:
            h = lab[v]
            if h == alpha[e]:
                if rng.random() < .05: Z.add((e, v))      # own coefficient may vanish
                else: keep.append(v)
            elif meas[h] and rng.random() < .9:           # nonzero restriction in R_h with c_o != 0
                clique.add(h); keep.append(v)
            else:
                Z.add((e, v))
                if meas[h] and rng.random() < .5: clique.add(h)   # nonzero restriction but zero own coeff
        for v in vs:                                       # extra foreign restrictions
            for h in range(g):
                if meas[h] and h != alpha[e] and rng.random() < .05: clique.add(h)
        for x in clique:
            for y in clique:
                if x != y: adj[x].add(y)
    # sanity: every edge of I_t has a measured end
    assert all(meas[x] or meas[y] for x in adj for y in adj[x])
    padj = {v: set() for v in range(n)}
    for e in cl:
        for x in e:
            for y in e:
                if x != y: padj[x].add(y)
    tau = tw(n, padj)
    M = [h for h in range(g) if meas[h]]; idx = {h: i for i, h in enumerate(M)}
    tadj = {i: set() for i in range(len(M))}
    for h in M:
        for y in adj[h]:
            if meas[y]: tadj[idx[h]].add(idx[y])
    for h in range(g):
        if not meas[h]:
            for x in adj[h]:
                for y in adj[h]:
                    if x != y: tadj[idx[x]].add(idx[y])
    Wa = wtw(len(M), tadj, [a[h] for h in M]); s_ = max(a)
    zv = len({v for (_, v) in Z})
    assert tau <= Wa + s_ - 1 + len(Z), (tau, Wa, s_, len(Z))
    cnt += 1; minslack = min(minslack, Wa + s_ - 1 + len(Z) - tau)
    mor_fail += not (Wa >= tau + 2 - s_ - len(Z))
print(f"Theorem C main bound: {cnt} random points pass, min slack {minslack}; 'moreover' (+2) fails in {mor_fail} random cases")
