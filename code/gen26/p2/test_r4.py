# Own exact check of R4: block enumeration (outer set S) + sweep, and divide & conquer, on K-cliques with lam=1.
import itertools, random, math, sys
random.seed(2)
def E(u, v): return (u, v) if u < v else (v, u)
def make_tables(K, lam=1):
    # W_w: dict from tuple of bits on edges at w (ordered by other endpoint) -> small signed int
    T = {}
    for w in range(K):
        nb = [x for x in range(K) if x != w]
        T[w] = {bits: random.randint(-3, 3) for bits in itertools.product(range(2), repeat=len(nb))}
    return T
def Wval(T, K, w, asg):
    return T[w][tuple(asg[E(w, x)] for x in range(K) if x != w)]
def brute(T, K):
    edges = [E(u, v) for u in range(K) for v in range(u+1, K)]; tot = 0
    for bits in itertools.product(range(2), repeat=len(edges)):
        asg = dict(zip(edges, bits)); p = 1
        for w in range(K): p *= Wval(T, K, w, asg)
        tot += p
    return tot
def sweep(T, K, fixed, order=None):
    """sweep in node order with edges in `fixed` (dict) fixed. Returns value, max #free strings enumerated in a step,
       max #free cut edges stored (table key length)."""
    order = order or list(range(K)); pos = {w: i for i, w in enumerate(order)}
    F = {(): 1}; cutkeys = []      # F maps assignment of free cut edges (tuple in order of cutkeys) -> value
    maxstep = 0; maxmem = 0; ops = 0
    for j in range(K):
        w = order[j]
        D = [e for e in cutkeys if w in e]                  # cut edges into w (to be summed out)
        Nn = [E(w, x) for x in order[j+1:] if E(w, x) not in fixed]   # new free edges
        newkeys = [e for e in cutkeys if w not in e] + Nn
        maxstep = max(maxstep, len(cutkeys) + len(Nn)); maxmem = max(maxmem, len(newkeys))
        G = {}
        for key, val in F.items():
            asg = dict(zip(cutkeys, key)); asg.update(fixed)
            for nb in itertools.product(range(2), repeat=len(Nn)):
                a2 = dict(asg); a2.update(zip(Nn, nb)); ops += 1
                # all edges at w are now assigned (earlier ones in cut or fixed, later in Nn or fixed)
                # edges at w to earlier nodes not in cut and not fixed cannot exist
                v = val*Wval(T, K, w, a2)
                k2 = tuple(a2[e] for e in newkeys)
                G[k2] = G.get(k2, 0) + v
        F = G; cutkeys = newkeys
    return F[()], maxstep, maxmem, ops
def outer_sweep(T, K, S):
    S = sorted(S); tot = 0; ops = 0; ms = mm = 0
    for bits in itertools.product(range(2), repeat=len(S)):
        v, ms, mm, o = sweep(T, K, dict(zip(S, bits))); tot += v; ops += o
    return tot, len(S)+ms, mm, ops
def dc(T, K, nodes, fixed, cnt):
    """value of sub-clique on `nodes` with outside edges fixed; D&C on middle cut. cnt[0] counts leaf terms."""
    k = len(nodes)
    if k <= 2:
        inner = [E(nodes[0], nodes[1])] if k == 2 else []
        tot = 0
        for bits in itertools.product(range(2), repeat=len(inner)):
            a = dict(fixed); a.update(zip(inner, bits)); p = 1; cnt[0] += 1
            for w in nodes: p *= Wval(T, K, w, a)
            tot += p
        return tot
    L, R = nodes[:k//2], nodes[k//2:]
    cut = [E(u, v) for u in L for v in R]; tot = 0
    for bits in itertools.product(range(2), repeat=len(cut)):
        f2 = dict(fixed); f2.update(zip(cut, bits))
        tot += dc(T, K, L, f2, cnt) * dc(T, K, R, f2, cnt)
    return tot
def dK(k): return 0 if k <= 1 else (k//2)*(k-k//2) + max(dK(k//2), dK(k-k//2))
for K in (4, 5, 6):
    for trial in range(3):
        T = make_tables(K); b = brute(T, K)
        v, ms, mm, ops = sweep(T, K, {})
        cnt = [0]; vd = dc(T, K, list(range(K)), {}, cnt)
        m = K//2; res = [f"K={K} brute={b} sweep={v} (step {ms}=s(K) {((K+1)**2)//4-1}, mem {mm}) dc={vd} (leaf terms {cnt[0]}, 2^d(K)={2**dK(K)})"]
        for d in range(1, m+1):
            S = [E(u, w) for u in range(m-d, m) for w in range(m, min(K, m+d))]
            vo, tm, mem, o = outer_sweep(T, K, S)
            res.append(f"   d={d}: outer+sweep={vo} time {tm} mem {mem}")
        print("\n".join(res) if trial == 0 else f"K={K} trial {trial}: brute={b} sweep={v} dc={vd} " +
              " ".join(f"d={d}:{outer_sweep(T, K, [E(u, w) for u in range(m-d, m) for w in range(m, min(K, m+d))])[0]}" for d in range(1, m+1)))
        assert b == v == vd
sys.stdout.flush()
