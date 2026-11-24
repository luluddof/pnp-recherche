"""Randomised balanced greedy (max-per-direction objective) for maximal consistent good sets.
Generic tau mod p, incremental elimination; best result re-verified exactly over Q with tau_k = k+1
and with random rational tau.  usage: python3 greedy_max.py l d trials seed"""
import sys, random, numpy as np
from core import edges, feasible_Q
from fractions import Fraction
p = 2**31 - 1
class Ech:
    def __init__(s, n): s.n = n; s.piv = []; s.rows = []
    def add(s, r):
        r = r.copy()
        for c, br in zip(s.piv, s.rows):
            if r[c]: r = (r - r[c] * br) % p
        nz = np.flatnonzero(r[:-1])
        if len(nz) == 0: return r[-1] == 0, False
        c = nz[0]; r = r * pow(int(r[c]), p - 2, p) % p
        s.piv.append(c); s.rows.append(r); return True, True
def run(l, d, trials, seed, strategy="rr"):
    rng = random.Random(seed); E = edges(l); n = (1 << l) * d
    best = None
    for tr in range(trials):
        taus = [rng.randrange(2, p - 1) for _ in range(l)]
        rowsE = []
        for (z, w, k) in E:
            r = np.zeros(n + 1, dtype=np.int64); pw = taus[k]
            for j in range(d):
                r[z * d + j] = (r[z * d + j] + pw) % p; r[w * d + j] = (r[w * d + j] + pw) % p; pw = pw * taus[k] % p
            r[n] = p - 2; rowsE.append(r)
        bydir = [[e for e, x in enumerate(E) if x[2] == k] for k in range(l)]
        for b in bydir: rng.shuffle(b)
        if strategy == "rr":
            order = []
            dirs = list(range(l))
            for i in range(len(bydir[0])):
                rng.shuffle(dirs); order += [bydir[k][i] for k in dirs]
        else:
            order = list(range(len(E))); rng.shuffle(order)
        ech = Ech(n); bad = []
        for e in order:
            ok, _ = ech.add(rowsE[e])
            if not ok: bad.append(e)
        per = [0] * l
        for e in bad: per[E[e][2]] += 1
        key = (max(per), len(bad))
        if best is None or key < best[0]: best = (key, per, bad)
    return best, E
if __name__ == "__main__":
    l, d, trials, seed = map(int, sys.argv[1:5]); strat = sys.argv[5] if len(sys.argv) > 5 else "rr"
    (mx, tot), per, bad, = run(l, d, trials, seed, strat)[0]; E = edges(l)
    good = [E[e] for e in range(len(E)) if e not in set(bad)]
    okQ = feasible_Q(l, d, good, [Fraction(k + 1) for k in range(l)]) if l <= 6 else None
    rr = random.Random(5); okQ2 = feasible_Q(l, d, good, [Fraction(rr.randrange(1, 1000), rr.randrange(1, 1000)) for _ in range(l)]) if l <= 6 else None
    print(f"l={l} d={d} {strat}: best max={mx} (rho_max={mx/2**l:.4f}) total={tot} per={per} exactQ(tau=k)={okQ} exactQ(rand)={okQ2}", flush=True)
