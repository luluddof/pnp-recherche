"""G4d experiment: IPS (Hilbert-like) certificate for an UNSAT 3-CNF built from
hypercube partial sums (telescoping multilinearisation), evaluated as an
oblivious arithmetic circuit (tensor-network contraction).

Certificate (over Q, or F_p with p > max occurrence):
  1 = sum_j c_j(y) prod_{i<j} (1 - c_i(y))            [clause part, size O(m)]
    + sum_k (y_k^2 - y_k) D_k(y)                         [Boolean part]
where c_j = falsifying indicator of clause j, P = prod_j (1 - c_j),
Q_k(y) = sum_{x_1..x_k in {0,1}} P(x_1..x_k, y_{k+1..n}) prod_{i<=k} l(y_i, x_i),
l(y,b) = y b + (1-y)(1-b), Q_0 = P, Q_n = ml(P) = 0 iff F UNSAT,
D_k = (Q_{k-1} - Q_k)/(y_k^2 - y_k)  (a polynomial: Q_k is Q_{k-1} linearly
interpolated in y_k). As a division-free circuit, D_k is obtained from
Q_{k-1} at d_k+1 integer values of y_k (d_k = #occurrences) by interpolation,
so size(D_k) <= (d_k+2) * size(Q_{k-1}) + O(d_k^2).

We check the identity at a random point mod p (= PIT, how IPS proofs are
verified), and count gates of Q_k computed (a) by brute force over 2^k
assignments, (b) by variable elimination along a min-degree order (cost
~ 2^treewidth). Claude's algorithm would replace (b) by 2^{O(m/sqrt log m)}
(asymptotically only).
"""
import random, sys, itertools

P = (1 << 61) - 1


def rand3cnf(n, m, rng):
    cls = set()
    while len(cls) < m:
        vs = rng.sample(range(n), 3)
        cls.add(tuple(sorted((v, rng.random() < 0.5) for v in vs)))
    return [list(c) for c in cls]


def satisfiable(n, cnf):
    for a in range(1 << n):
        if all(any(((a >> v) & 1) == s for v, s in c) for c in cnf):
            return True
    return False


def lit(val, sign):  # value of literal (x if sign else 1-x), val in F_p
    return val % P if sign else (1 - val) % P


def clause_factor(c, summed, y):
    """table over the summed variables of clause c of 1 - c_j."""
    sv = [v for v, _ in c if v in summed]
    tab = {}
    for bits in itertools.product((0, 1), repeat=len(sv)):
        asg = dict(zip(sv, bits))
        f = 1
        for v, s in c:
            val = asg[v] if v in summed else y[v]
            f = f * (1 - lit(val, s)) % P
        tab[bits] = (1 - f) % P
    return (tuple(sv), tab)


def Q_value(cnf, n, order_k, y, method):
    """Q over summed set = order_k (list of vars) at point y; returns (value, gates)."""
    summed = set(order_k)
    factors = [clause_factor(c, summed, y) for c in cnf]
    for v in order_k:  # l(y_v, x_v)
        factors.append(((v,), {(0,): (1 - y[v]) % P, (1,): y[v] % P}))
    gates = 0
    if method == "brute":
        k = len(order_k)
        tot = 0
        for bits in itertools.product((0, 1), repeat=k):
            asg = dict(zip(order_k, bits))
            prod = 1
            for sc, tab in factors:
                prod = prod * tab[tuple(asg[v] for v in sc)] % P
            gates += len(factors)
            tot = (tot + prod) % P
            gates += 1
        return tot, gates
    # variable elimination, min-degree order
    remaining = set(order_k)
    while remaining:
        def deg(v):
            s = set()
            for sc, _ in factors:
                if v in sc:
                    s |= set(sc)
            return len(s)
        v = min(remaining, key=deg)
        remaining.discard(v)
        inv = [f for f in factors if v in f[0]]
        factors = [f for f in factors if v not in f[0]]
        scope = sorted(set().union(*[set(sc) for sc, _ in inv]) - {v})
        tab = {}
        for bits in itertools.product((0, 1), repeat=len(scope)):
            asg = dict(zip(scope, bits))
            tot = 0
            for b in (0, 1):
                asg[v] = b
                prod = 1
                for sc, t in inv:
                    prod = prod * t[tuple(asg[u] for u in sc)] % P
                gates += len(inv)
                tot = (tot + prod) % P
                gates += 1
            tab[bits] = tot
        factors.append((tuple(scope), tab))
    val = 1
    for sc, t in factors:
        val = val * t[()] % P
        gates += 1
    return val, gates


def run(n, ratio, seed, methods=("ve",)):
    rng = random.Random(seed)
    m = int(round(ratio * n))
    while True:
        cnf = rand3cnf(n, m, rng)
        if n > 22 or not satisfiable(n, cnf):
            break
    y = [rng.randrange(2, P) for _ in range(n)]
    order = list(range(n))
    occ = [sum(1 for c in cnf for v, _ in c if v == i) for i in range(n)]
    out = {}
    for meth in methods:
        Qs, gates = [], 0
        for k in range(n + 1):
            val, g = Q_value(cnf, n, order[:k], y, meth)
            Qs.append(val)
            if k >= 1:
                gates += (occ[order[k - 1]] + 2) * g  # interpolation points
        # clause part
        cs = []
        for c in cnf:
            f = 1
            for v, s in c:
                f = f * (1 - lit(y[v], s)) % P
            cs.append(f)
        clause_part, pref = 0, 1
        for cj in cs:
            clause_part = (clause_part + cj * pref) % P
            pref = pref * (1 - cj) % P
        bool_part = 0
        for k in range(1, n + 1):
            v = order[k - 1]
            q = (y[v] * y[v] - y[v]) % P
            Dk = (Qs[k - 1] - Qs[k]) * pow(q, P - 2, P) % P
            bool_part = (bool_part + q * Dk) % P
        ok = (clause_part + bool_part) % P == 1
        out[meth] = dict(Qn_zero=Qs[n] == 0, identity_holds=ok, gates=gates)
    return m, out


if __name__ == "__main__":
    import math
    ratio = 4.26
    print("n m method Q_n==0 identity log2(cert_gates)")
    for n in [8, 10, 12, 14, 16, 18, 20, 24, 28, 32, 40, 48]:
        meths = ("brute", "ve") if n <= 16 else ("ve",)
        m, out = run(n, ratio, seed=n, methods=meths)
        for meth, r in out.items():
            print(n, m, meth, r["Qn_zero"], r["identity_holds"],
                  round(math.log2(r["gates"]), 2), flush=True)
