#!/usr/bin/env python3
"""G18e: numerical search for DF decompositions on tiny instances (r = 2 groups per direction).
instance: parties = ints (labels in {0,1}^u), wires = list of (p,q); direction i groups parties by bit i.
Find F_k, G_k with sum_k F_k(bits in group 0 of dir d_k) G_k(bits in group 1) = I.  Real numerics (LM),
this is EVIDENCE only; exact decompositions are verified separately in exact.py."""
import sys, itertools, numpy as np
from scipy.optimize import least_squares

def setup(wires, u):
    M = len(wires); nb = 2 * M; T = 1 << nb
    party = []
    for (p, q) in wires: party += [p, q]
    idx = np.arange(T); bits = [(idx >> b) & 1 for b in range(nb)]
    I = np.ones(T)
    for e in range(M): I *= (bits[2 * e] == bits[2 * e + 1])
    dirs = []
    for i in range(u):
        g0 = [b for b in range(nb) if (party[b] >> i) & 1 == 0]; g1 = [b for b in range(nb) if (party[b] >> i) & 1 == 1]
        a = sum(bits[b] << k for k, b in enumerate(g0)) if g0 else np.zeros(T, dtype=int)
        c = sum(bits[b] << k for k, b in enumerate(g1)) if g1 else np.zeros(T, dtype=int)
        dirs.append((a, 1 << len(g0), c, 1 << len(g1)))
    return I, dirs

def solve(wires, u, prof, tries=20, seed=0, tol=1e-18, maxnfev=400):
    I, dirs = setup(wires, u); T = len(I)
    terms = [i for i in range(u) for _ in range(prof[i])]
    offs = []; n = 0
    for d in terms:
        a, na, c, nc = dirs[d]; offs.append((n, n + na)); n += na + nc
    def unpack(x):
        return [(x[o1:o2], x[o2:o2 + dirs[d][3]]) for (o1, o2), d in zip(offs, terms)]
    def res(x):
        r = -I.copy()
        for (F, G), d in zip(unpack(x), terms):
            a, na, c, nc = dirs[d]; r += F[a] * G[c]
        return r
    rows = np.arange(T)
    def jac(x):
        J = np.zeros((T, n))
        for (F, G), d, (o1, o2) in zip(unpack(x), terms, offs):
            a, na, c, nc = dirs[d]
            J[rows, o1 + a] = G[c]; J[rows, o2 + c] = F[a]
        return J
    rng = np.random.default_rng(seed); best = (np.inf, None)
    for t in range(tries):
        x0 = rng.normal(size=n) * 0.7
        s = least_squares(res, x0, jac=jac, method='trf', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=maxnfev)
        c = 2 * s.cost
        if c < best[0]: best = (c, s.x)
        if c < tol: break
    return best[0], best[1], t + 1

if __name__ == "__main__":
    case = sys.argv[1]
    if case == "pk2":      # one full packet a=u=2: slots 00,01,10,11 ; wire of slot x joins x and ~x
        wires = [(0, 3), (3, 0), (1, 2), (2, 1)]; u = 2; profs = [(2,2),(1,3),(2,3),(1,4),(3,3),(2,4),(1,5),(0,7),(3,4),(2,5),(1,6)]
    elif case == "pk2m3":
        wires = [(0, 3), (3, 0), (1, 2)]; u = 2; profs = [(2, 1), (1, 2), (3, 0), (2, 2)]
    elif case == "pk3m3":  # a=u=3, three slots on three different antipodal pairs
        wires = [(0, 7), (1, 6), (2, 5)]; u = 3; profs = [p for N in (2, 3) for p in itertools.product(range(N + 1), repeat=3) if sum(p) == N] + [(2,2,0),(2,0,2),(0,2,2)]
    elif case == "pk3m4":  # a=u=3, one slot on each of the four antipodal pairs
        wires = [(0, 7), (1, 6), (2, 5), (3, 4)]; u = 3
        profs = [p for N in (2, 3, 4) for p in itertools.product(range(N + 1), repeat=3) if sum(p) == N and p[0] >= p[1] >= p[2]]
        profs += [p for N in (5, 6, 7) for p in itertools.product(range(N + 1), repeat=3) if sum(p) == N and p[0] >= p[1] >= p[2]]
    elif case == "k4":     # labels {0,1}^2, one wire per unordered pair of distinct labels
        wires = list(itertools.combinations(range(4), 2)); u = 2
        profs = [(N // 2, N - N // 2) for N in range(4, 16)] + [(k, N - k) for N in (8, 10, 12) for k in (0, 1, 2, 3)]
    tries = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    for prof in profs:
        c, x, t = solve(wires, u, prof, tries=tries)
        print(f"{case} profile {prof} N={sum(prof)}: best residual^2 = {c:.3e} after {t} starts", "FOUND" if c < 1e-18 else "", flush=True)
