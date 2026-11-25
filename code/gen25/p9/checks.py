# G25-P9 small exact checks.
# (1) tw(L(K_n)) exactly (subset DP over elimination orders) vs s(n)-1, s(n) = floor((n+1)^2/4)-1 (sweep bag size).
# (2) Row sweep on a general network (Definition 2 type: node settings z_a, edge factors that factor through
#     c_ab bits) equals brute force; tiny random instances.
# (3) Random-order cut bound of Lemma RS: max_j cut_j of the best of a few random orders vs U/4 and the bound.
# (4) r~(lam) := max_s sqrt(c_row^2 + c_col^2)/(2 lam) for App. B (HR leading term, as p2's rho_bar), vs sqrt2*rho_bar.
import itertools, random, math, sys
from functools import lru_cache
sys.path.insert(0, '.')

def tw_exact(n_vert, adj):
    full = (1 << n_vert) - 1
    def Q(S, v):   # vertices outside S+{v} reachable from v through S
        seen = 1 << v; stack = [v]; out = 0
        while stack:
            u = stack.pop()
            nb = adj[u] & ~seen
            seen |= nb
            m = nb
            while m:
                b = m & -m; w = b.bit_length() - 1; m ^= b
                if (S >> w) & 1: stack.append(w)
                else: out |= b
        return bin(out).count('1')
    TW = {0: -1}
    for size in range(1, n_vert + 1):
        for comb in itertools.combinations(range(n_vert), size):
            S = 0
            for c in comb: S |= 1 << c
            best = 10**9
            for v in comb:
                Sv = S & ~(1 << v)
                val = max(TW[Sv], Q(Sv, v))
                if val < best: best = val
            TW[S] = best
    return TW[full]

def line_graph_complete(n):
    E = [(a, b) for a in range(n) for b in range(a+1, n)]
    adj = [0]*len(E)
    for i, e in enumerate(E):
        for j, f in enumerate(E):
            if i != j and set(e) & set(f): adj[i] |= 1 << j
    return len(E), adj

print("(1) treewidth of the line graph of K_n (exact) vs sweep bag size s(n) = floor((n+1)^2/4)-1")
for n in range(3, 7):
    nv, adj = line_graph_complete(n)
    tw = tw_exact(nv, adj); s = ((n+1)**2)//4 - 1
    print(f"   n={n}: tw(L(K_n)) = {tw}, tw+1 = {tw+1}, s(n) = {s}")
sys.stdout.flush()

# (2) general network sweep vs brute force
def check_network(rng, g, cmax, P):
    # nodes 0..g-1, settings z_a in [0, 2^P), edge {a,b}: theta in [0, 2^{c_ab}), factors f_ab(z_a, theta), f_ba(z_b, theta)
    c = {(a, b): rng.randint(0, cmax) for a in range(g) for b in range(a+1, g)}
    W = [[rng.randint(-2, 2) for _ in range(2**P)] for a in range(g)]
    f = {}
    for (a, b), cab in c.items():
        f[(a, b)] = [[rng.randint(-1, 1) for _ in range(2**cab)] for _ in range(2**P)]
        f[(b, a)] = [[rng.randint(-1, 1) for _ in range(2**cab)] for _ in range(2**P)]
    # brute force: value = sum_z prod W_a(z_a) prod_{a<b} Psi_ab(z_a,z_b), Psi_ab = sum_theta f_ab f_ba
    tot = 0
    for z in itertools.product(range(2**P), repeat=g):
        pr = 1
        for a in range(g): pr *= W[a][z[a]]
        if pr == 0: continue
        for (a, b), cab in c.items():
            pr *= sum(f[(a, b)][z[a]][th]*f[(b, a)][z[b]][th] for th in range(2**cab))
            if pr == 0: break
        tot += pr
    # sweep in a random order: table over theta's on the cut edges; node table Wt_a(thetas) = sum_z W f...
    order = list(range(g)); rng.shuffle(order); pos = {v: i for i, v in enumerate(order)}
    def key(a, b): return (a, b) if a < b else (b, a)
    F = {(): 1}; cut = []
    for j, v in enumerate(order):
        inc = [e for e in cut if v in e]; keep = [e for e in cut if v not in e]
        new = [key(v, u) for u in order[j+1:]]
        G = {}
        for k, val in F.items():
            s = dict(zip(cut, k))
            for nv in itertools.product(*[range(2**c[e]) for e in new]):
                sn = dict(zip(new, nv))
                # node v's table at these edge strings
                wt = 0
                for zv in range(2**P):
                    t = W[v][zv]
                    if t == 0: continue
                    for e in inc:
                        other = e[0] if e[1] == v else e[1]
                        t *= f[(v, other)][zv][s[e]]
                    for e in new:
                        other = e[0] if e[1] == v else e[1]
                        t *= f[(v, other)][zv][sn[e]]
                    wt += t
                if wt:
                    kk = tuple(s[e] for e in keep) + nv
                    G[kk] = G.get(kk, 0) + val*wt
        F = G; cut = keep + new
    return tot, F.get((), 0)

rng = random.Random(7); ok = True
for trial in range(12):
    g = rng.choice([3, 4, 5]); P = 1 if g == 5 else rng.choice([1, 2]); cmax = 2 if g < 5 else 1
    b, s = check_network(rng, g, cmax, P)
    ok &= (b == s)
print(f"(2) row/network sweep == brute force on 12 random Definition-2 networks: {'ALL OK' if ok else 'MISMATCH'}")
sys.stdout.flush()

# (3) cut bound: weights c_ab in [0, Pg] (Pg = P/g), adversarial-ish structures
import numpy as np
def maxcut_of_order(Cm, order):
    g = len(order); inS = np.zeros(g, bool); cur = 0.0; mx = 0.0
    for v in order:
        # adding v: cut changes by (weights v->outside) - (weights v->inside)
        out = Cm[v][~inS].sum() - Cm[v][v]; ins = Cm[v][inS].sum()
        cur += out - ins; inS[v] = True; mx = max(mx, cur)
    return mx
nrng = np.random.default_rng(3)
print("(3) Lemma RS: best of 20 random orders, max_j cut_j / (U/4) and the proved bound / (U/4)")
for g in (64, 256, 1024):
    for kind in ('uniform', 'random', 'two-block', 'heavy-half'):
        if kind == 'uniform': Cm = np.ones((g, g))
        elif kind == 'random': Cm = nrng.random((g, g)); Cm = (Cm + Cm.T)/2
        elif kind == 'two-block':
            Cm = np.zeros((g, g)); h = g//2; Cm[:h, h:] = 1; Cm[h:, :h] = 1
        else:
            Cm = np.full((g, g), 0.2); Cm[:g//2, :g//2] = 1
        np.fill_diagonal(Cm, 0)
        U = Cm.sum(); Pmax = Cm.sum(1).max(); Pg = Cm.max()
        best = min(maxcut_of_order(Cm, nrng.permutation(g)) for _ in range(20))
        # bound: U/4 * g/(g-1) + Pg*g * sqrt(ln(2g^2+2g+1) g/4) (1+1/g)  [P = g*Pg]
        bound = U/4*g/(g-1) + g*Pg*math.sqrt(math.log(2*g*g+2*g+1)*g/4)*(1+1/g)
        print(f"   g={g:5d} {kind:10s}: best/(U/4) = {best/(U/4):.4f}   bound/(U/4) = {bound/(U/4):.4f}")
sys.stdout.flush()

# (4) r~(lam)
import thmA_explicit as T
print("(4) App. B frontier: sqrt2*rho_bar(lam) [= max_s (c_row+c_col)/(2 sqrt2 lam)] and r~(lam) = max_s sqrt(c_row^2+c_col^2)/(2 lam)")
for L2 in (10, 12, 16, 20, 30, 40):
    lam = 2**L2; h = T.h_of(lam)
    best_r = 0; best_sum = 0
    for s in range(-2*h, 2*h + 1, max(1, (4*h)//4000)):
        cr = min(2*lam, math.ceil(math.log2(2*h+1) + T.lp_asym(max(2*h+s, 1)) + 0.5))
        cc = min(2*lam, math.ceil(math.log2(2*h+1) + T.lp_asym(max(2*h-s, 1)) + 0.5))
        best_r = max(best_r, math.hypot(cr, cc)/(2*lam)); best_sum = max(best_sum, (cr+cc)/(2*math.sqrt(2)*lam))
    print(f"   lam=2^{L2}: sqrt2*rho_bar = {best_sum:.5f}   r~ = {best_r:.5f}")
