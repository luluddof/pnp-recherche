# Scan: lattice VOAs (rank 1-4), modules lam+L, regradings xi in the Voronoi cell; per-copy (tensor-power family) criterion.
#  phi^{(z)}_X(b) = log sum_{p in X} exp(-b(|p|^2/2 - <z,p>)) - r sum_k log(1-e^{-bk});  s(x) = inf_b [phi(b) + b x]
#  columns:  LOOSE  = labels & pairs both in the cost grading (G22e "loose"; upper bound on N for every construction)
#            B0     = branch (i)/channel 0: pairs in the CONFORMAL grading of V, labels in M^xi
#            SC     = best single channel g in S_xi (and two-channel mixtures 0/g, pi in {1/4,1/2,3/4}), rigorous Conf-type bounds
import itertools, math, sys, numpy as np
SQ8 = 2*math.sqrt(2)
bs = np.exp(np.linspace(math.log(25), math.log(0.08), 500))
def lattice_norms(G, coset, zeta, R2):
    G = np.array(G, float); r = len(G)
    ev = np.linalg.eigvalsh(G).min(); B = int(math.ceil(math.sqrt(R2/ev))) + 2
    rng = np.arange(-B, B+1)
    pts = np.array(list(itertools.product(rng, repeat=r)), float) + np.array(coset, float)
    nrm = np.einsum('ij,jk,ik->i', pts, G, pts)
    pts, nrm = pts[nrm <= R2], nrm[nrm <= R2]
    w = nrm/2 - pts @ G @ np.array(zeta, float)
    return w
def phi(G, coset, zeta, R2):
    r = len(G); w = lattice_norms(G, coset, zeta, R2)
    k = np.arange(1, 4000)
    out = []
    for b in bs:
        t = -b*w; m = t.max()
        out.append(m + math.log(np.exp(t-m).sum()) - r*np.log1p(-np.exp(-b*k)).sum())
    return np.array(out)
def rate(P, x): return float(np.min(P + bs*x))
def crit(Pcost, lam_fns, grid):
    best = (9, 0, 0)
    for u in grid:
        for v in grid:
            lam = min(f(u, v) for f in lam_fns)
            if lam > 1e-6:
                c = 2*rate(Pcost, u+v)/lam
                if c < best[0]: best = (c, u, v)
    return best
def S_xi(G, xi):
    G = np.array(G, float); r = len(G)
    out = []
    for n in itertools.product(range(-3, 4), repeat=r):
        g = np.array(n, float)
        if abs(g@G@g/2 - g@G@np.array(xi)) < 1e-9: out.append(tuple(n))
    return out
def run(name, G, xi, coset=None, R2=None, rank_scale=1.0):
    r = len(G); coset = coset or [0]*r
    R2 = R2 or (60 if r <= 2 else 36 if r == 3 else 22)
    xi = np.array(xi, float)
    S = S_xi(G, xi)
    if lattice_norms(G, [0]*r, xi, 30).min() < -1e-9: print("%-34s xi=%s OUTSIDE Voronoi cell, skipped" % (name, np.round(xi,3).tolist())); return None
    PC = phi(G, coset, xi, R2)               # cost profile of M (grading xi)
    PVx = phi(G, [0]*r, xi, R2)              # V in grading xi (loose pairs)
    PV0 = phi(G, [0]*r, [0]*r, R2)           # V conformal (channel 0 pairs)
    grid = np.exp(np.linspace(math.log(0.03*r), math.log(3.0*r), 38))
    U = lambda u, v: rate(PC, 2*u+2*v)/2
    lab = lambda u, v: rate(PC, u)
    L = crit(PC, [lab, lambda u, v: rate(PVx, v), U], grid)
    B0 = crit(PC, [lab, lambda u, v: rate(PV0, v), U], grid)
    sc = (9, None)
    for g in S:
        if all(x == 0 for x in g): continue
        g = np.array(g, float)
        PVg = phi(G, [0]*r, g/2, R2); PMg = phi(G, coset, xi - g/2, R2)
        PM0 = PC
        for pi in [0.25, 0.5, 0.75, 1.0]:
            Pr = (1-pi)*PV0 + pi*PVg; Pl = (1-pi)*PM0 + pi*PMg
            c = crit(PC, [lab, lambda u, v: rate(Pr, v), lambda u, v: rate(Pl, u), U], grid)
            if c[0] < sc[0]: sc = (c[0], (tuple(int(x) for x in g), pi), c[1], c[2])
    tag = "DIP-LOOSE<2sqrt2" if L[0] < SQ8 else ""
    print("%-34s xi=%-22s |S|=%d  LOOSE %.4f (e=%.2f,v=%.2f)  B0 %.4f  SC %s  %s" % (
        name, np.round(xi, 3).tolist(), len(S), L[0], L[1], L[2], B0[0],
        ("%.4f %s" % (sc[0], sc[1])) if sc[1] else "  -   ", tag), flush=True)
    return L, B0, sc
if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    A1 = [[2]]; A2 = [[2, -1], [-1, 2]]; A1A1 = [[2, 0], [0, 2]]
    A3 = [[2, -1, 0], [-1, 2, -1], [0, -1, 2]]; D4 = [[2, -1, 0, 0], [-1, 2, -1, -1], [0, -1, 2, 0], [0, -1, 0, 2]]
    if which in ("all", "r1"):
        for t in [0, 0.25, 0.5, 0.75, 0.9, 1.0]:
            run("A1 V", A1, [t/2])
        for t in [0, 0.5, 1.0]:
            run("A1 M=V_{a/2+L}", A1, [t/2], coset=[0.5])
        run("A1(x)2 = sqrt(4)Z? no: Z (odd, VOSA)", [[1]], [0.5])
    if which in ("all", "r2"):
        # A2: Voronoi cell vertices = fundamental weights (deep holes), edge midpoints = roots/2
        for xi in [[0, 0], [0.5, 0], [1/3, 2/3], [2/3, 1/3], [0.25, 0.1], [0.4, 0.3], [0.5, 0.5]]:
            run("A2 V", A2, xi)
        for xi in [[0, 0], [1/3, 2/3], [0.5, 0]]:
            run("A2 M=V_{w1+L}", A2, xi, coset=[2/3, 1/3])
        for xi in [[0.5, 0], [0.5, 0.5], [0.25, 0.5]]:
            run("A1+A1 V", A1A1, xi)
    if which in ("all", "r3"):
        for xi in [[0.5, 0, 0], [0.5, 0.5, 0.5], [0.5, 1, 0.5], [0.25, 0.5, 0.75]]:
            run("A3 V", A3, xi)
    if which in ("all", "r4"):
        for xi in [[0.5, 0, 0, 0], [0.5, 1, 0.5, 0.5], [1, 1.5, 1, 1]]:
            run("D4 V", D4, xi)
