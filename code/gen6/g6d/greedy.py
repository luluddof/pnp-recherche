"""Randomised greedy maximal-consistent good-edge sets (commutative class), mod-p incremental
elimination, final exact check over Q (milp.exact_feasible).
Reformulation used: r_z in P_{d-1}; edge (z,z^k) good iff (r_z + r_{z^k})(tau_k) = 2/tau_k.
(Equivalent to g_z = chi(z) - t chi(z) r_z(t), deg g <= d, g_z(0)=chi(z).)
usage: python3 greedy.py l d trials [order: random|layer]"""
import sys, itertools, random
import numpy as np
P = 2147483629  # prime < 2^31

def inv(a): return pow(int(a), P - 2, P)

class Echelon:
    def __init__(self, n):
        self.n = n; self.piv = {}  # pivot col -> row (np array length n+1)
    def add(self, row):
        row = row.copy() % P
        for c in np.nonzero(row[:-1])[0]:
            pass
        while True:
            nz = np.nonzero(row[:-1])[0]
            nz = [c for c in nz if c in self.piv]
            if not nz: break
            c = nz[0]
            row = (row - row[c] * self.piv[c]) % P
        nz = np.nonzero(row[:-1])[0]
        if len(nz) == 0:
            return row[-1] == 0          # consistent iff rhs 0
        c = nz[0]
        row = (row * inv(row[c])) % P
        # keep basis reduced w.r.t. new pivot
        for c2, r2 in self.piv.items():
            if r2[c]: self.piv[c2] = (r2 - r2[c] * row) % P
        self.piv[c] = row
        return True

def build(l, d):
    taus = list(range(1, l + 1))
    V = list(itertools.product([0, 1], repeat=l)); idx = {v: i for i, v in enumerate(V)}
    edges = []
    for w in V:
        for k in range(l):
            if w[k] == 0:
                u = list(w); u[k] = 1; edges.append((w, tuple(u), k))
    n = len(V) * d
    rows = []
    for (w, u, k) in edges:
        row = np.zeros(n + 1, dtype=np.int64)
        tk = taus[k]
        for j in range(d):  # r(t) = sum_j c_j t^j
            row[idx[w] * d + j] = (row[idx[w] * d + j] + pow(tk, j, P)) % P
            row[idx[u] * d + j] = (row[idx[u] * d + j] + pow(tk, j, P)) % P
        row[-1] = (2 * inv(tk)) % P
        rows.append(row)
    return V, edges, rows, n

def greedy(l, d, order, data):
    V, edges, rows, n = data
    E = Echelon(n); bad = []
    for e in order:
        if not E.add(rows[e]): bad.append(e)
    return bad

def run(l, d, trials, mode="random", seed=1):
    data = build(l, d); V, edges, rows, n = data
    rng = random.Random(seed)
    best = None
    for tr in range(trials):
        order = list(range(len(edges)))
        rng.shuffle(order)
        if mode == "layer":  # process low-weight lower endpoints first
            order.sort(key=lambda e: (sum(edges[e][0]) + rng.random() * 2))
        bad = greedy(l, d, order, data)
        per = [0] * l
        for e in bad: per[edges[e][2]] += 1
        key = (len(bad), max(per))
        if best is None or key < best[0]: best = (key, bad, per)
    return best, data

if __name__ == "__main__":
    l, d, trials = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    mode = sys.argv[4] if len(sys.argv) > 4 else "random"
    (key, bad, per), data = run(l, d, trials, mode)
    tot = len(bad)
    print(f"l={l} d={d} trials={trials} {mode}: best_total={tot} per_dir={per} frac={tot/(l*2**(l-1)):.4f} "
          f"max_rho={max(per)/2**l:.4f} conj={(l-d)*2**(l-d-1)}", flush=True)
