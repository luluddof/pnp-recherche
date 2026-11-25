# Own toy implementation of R1: derandomised random clustering (method of conditional expectations).
import numpy as np, math
rng = np.random.default_rng(3)
m = 6000; n = int(m/4.26); g = 6
cl = np.array([rng.choice(n, 3, replace=False) for _ in range(m)])
occ = {}
for c_, vs in enumerate(cl):
    for v in vs: occ.setdefault(v, []).append(c_)
edges = [(lst[i], lst[i+1]) for lst in occ.values() for i in range(len(lst)-1) if lst[i] != lst[i+1]]
W = len(edges)
# greedy proper edge colouring (<= 2*6-1 = 11 colours; Shannon's 9 not needed for the test)
used = [set() for _ in range(m)]; col = []
for u, v in edges:
    k = 0
    while k in used[u] or k in used[v]: k += 1
    used[u].add(k); used[v].add(k); col.append(k)
col = np.array(col); C = col.max()+1
print(f"m={m} n={n} W={W} colours={C} maxdeg={max(len(s) for s in used)}")
pairs = [(a, b) for a in range(g) for b in range(a+1, g)]
mu_bar = 2*W/g**2      # expected c_ab ... per unordered pair: W * 2/g^2
cnt = np.bincount(col, minlength=C)
def setup(d):
    # per event (pair, colour i): thr = mu_i + d (mu_i + mu_bar/C), s = ln(thr/mu_i)
    mu_i = cnt*2/g**2; thr = mu_i + d*(mu_i + mu_bar/C); s = np.log(thr/mu_i)
    # sizes: thr_c = (1+dc) m/g, s_c = ln(1+dc)
    return mu_i, thr, s
def logfactor(st_u, st_v, a, b, s):   # log E[exp(s 1[{lu,lv}={a,b}]) | partial labels], -1 = unlabelled
    es = math.exp(s)
    if st_u >= 0 and st_v >= 0: return s if {st_u, st_v} == {a, b} else 0.0
    if st_u >= 0 or st_v >= 0:
        l = st_u if st_u >= 0 else st_v
        return math.log(1 + (es-1)/g) if l in (a, b) else 0.0
    return math.log(1 + (es-1)*2/g**2)
def pairpart(d):
    mu_i, thr, s = setup(d); return len(pairs)*np.exp(-s*thr + cnt*np.log(1 + (np.exp(s)-1)*2/g**2)).sum()
def sizepart(dc):
    sc = math.log(1+dc); return g*math.exp(-sc*(1+dc)*m/g + m*math.log(1 + (math.exp(sc)-1)/g))
lo, hi = 0.0, 2.0
for _ in range(60):
    mid = (lo+hi)/2; lo, hi = (mid, hi) if pairpart(mid) > 0.25 else (lo, mid)
d = hi; lo, hi = 0.0, 2.0
for _ in range(60):
    mid = (lo+hi)/2; lo, hi = (mid, hi) if sizepart(mid) > 0.25 else (lo, mid)
dc = hi
mu_i, thr, s = setup(d); sc = math.log(1+dc); thc = (1+dc)*m/g
lab = -np.ones(m, int)
# log-terms per pair event (p, i) and size events
LT = np.zeros((len(pairs), C)); 
for p, (a, b) in enumerate(pairs):
    LT[p] = -s*thr + cnt*np.log(1 + (np.exp(s)-1)*2/g**2)
LS = np.full(g, -sc*thc + m*math.log(1 + (math.exp(sc)-1)/g))
Phi = lambda: np.exp(LT).sum() + np.exp(LS).sum()
print(f"d={d} dc={dc}: Phi(empty) = {Phi():.4f}")
inc = [[] for _ in range(m)]
for k, (u, v) in enumerate(edges): inc[u].append(k); inc[v].append(k)
hist = [Phi()]; mono_ok = True
for u in range(m):
    best = None
    for l in range(g):
        dLT = np.zeros_like(LT)
        for k in inc[u]:
            a_, b_ = edges[k]; o = b_ if a_ == u else a_; i = col[k]
            for p, (a, b) in enumerate(pairs):
                dLT[p, i] += logfactor(l, lab[o], a, b, s[i]) - logfactor(-1, lab[o], a, b, s[i])
        dLS = np.zeros(g); dLS[:] = -math.log(1 + (math.exp(sc)-1)/g); dLS[l] += sc
        val = np.exp(LT+dLT).sum() + np.exp(LS+dLS).sum()
        if best is None or val < best[0]: best = (val, l, dLT, dLS)
    val, l, dLT, dLS = best; LT += dLT; LS += dLS; lab[u] = l
    if val > hist[-1]*(1+1e-12): mono_ok = False
    hist.append(val)
# verify events
fail = 0; maxratio = 0
for i in range(C):
    sel = np.where(col == i)[0]
    for p, (a, b) in enumerate(pairs):
        X = sum(1 for k in sel if {lab[edges[k][0]], lab[edges[k][1]]} == {a, b})
        fail += X >= thr[i]
sizes = np.bincount(lab, minlength=g)
print(f"Phi final {hist[-1]:.4f}; monotone nonincreasing: {mono_ok}; failed pair-colour events {fail}; max size {sizes.max()} vs thr {thc:.1f}")
# overall loads
cab = np.zeros((g, g))
for (u, v) in edges:
    if lab[u] != lab[v]: cab[lab[u], lab[v]] += 1; cab[lab[v], lab[u]] += 1
print(f"max c_ab = {cab.max():.0f}, mu_bar = {mu_bar:.1f}, (1+2d) mu_bar = {(1+2*d)*mu_bar:.1f}")
