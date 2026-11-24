# Exact profiles (generating functions) for pairs (V, M) in the GENERALISED accounting of Appendix B (see notes.md sect. 1):
#   N = D_V(h) pairs of locality order <= H = 2h (h in the charge-symmetric grading of V),
#   labels = N basis vectors of M_{<=E0}, E0 least with D_M(E0) >= N,  budget B = 2 E0 + H,
#   cost(x) = L_M(x) + L_M(B - x), L = log2 D (cumulative dims);  ratio = max_x cost / (2 log2 N)   (circle: sqrt2; line: 1; rho = ratio/2)
# Weights are stored in integer units of 1/U.
import math, sys
def series(towers, W, pref=1):
    """towers: list of (kind, start, step, mult) in integer units; kind 'b' boson, 'f' fermion. Returns coefficients up to W."""
    c = [0]*(W+1); c[0] = pref
    for kind, start, step, mult in towers:
        r = start
        while r <= W:
            for _ in range(mult):
                if kind == 'f':
                    for w in range(W, r-1, -1): c[w] += c[w-r]
                else:
                    for w in range(r, W+1): c[w] += c[w-r]
            r += step
    return c
def cum(c):
    out = []; s = 0
    for v in c: s += v; out.append(s)
    return out
def analyse(name, DV, DM, U, hlist, verbose=False):
    best = None
    for h in hlist:            # h in units 1/U ; H = 2h
        N = DV[h]
        if N < 4: continue
        E0 = next((E for E in range(len(DM)) if DM[E] >= N), None)
        if E0 is None: continue
        B = 2*E0 + 2*h
        if B >= len(DM): break
        lam = math.log2(N)
        costs = [math.log2(DM[x]) + math.log2(DM[B-x]) for x in range(B+1)]
        r = max(costs)/(2*lam); rmin = min(costs)/(2*lam)
        rec = (r, h/U, E0/U, lam, rmin, math.log2(DM[0])/lam)
        if verbose: print(f"   h={h/U:7.2f} E0={E0/U:7.2f} lam={lam:9.2f} ratio={r:.4f} (min over x: {rmin:.4f}; L_M(0)/lam={rec[5]:.3f})")
        if best is None or r < best[0]: best = rec
    if best: print(f"{name:52s} best ratio {best[0]:.4f} at h={best[1]:.2f} E0={best[2]:.2f} lam={best[3]:.1f}  [line check min_x={best[4]:.3f}, L_M(0)/lam={best[5]:.3f}]", flush=True)
    return best
if __name__ == "__main__":
    U = 2
    Wmax = int(sys.argv[1]) if len(sys.argv) > 1 else 400   # in units 1/2
    hl = list(range(1, Wmax//4))
    print("ratio = max_x[L_M(x)+L_M(B-x)]/(2 lam); circle sqrt2 = 1.4142; concave-(V) bound 4/3; line 1")
    # 1. Heisenberg (App B)
    D = cum(series([('b', 2, 2, 1)], 4*Wmax)); analyse("Heisenberg rank 1 (App B), M=V", D, D, U, list(range(2, Wmax, 2)))
    for d in (2, 8, 24, 64):
        D = cum(series([('b', 2, 2, d)], Wmax)); analyse(f"Heisenberg rank {d}, M=V", D, D, U, hl)
        Dt = cum(series([('b', 1, 2, d)], Wmax)); analyse(f"Heisenberg rank {d}, M=Z2-twisted (half-odd modes)", D, Dt, U, hl)
        Dt2 = [v*2**(d//2) for v in Dt]; analyse(f"lattice-type twisted, ground 2^(d/2), V=Heis^{d}", D, Dt2, U, hl)
    for n in (8, 16, 32, 64, 128):
        NS = cum(series([('f', 1, 2, n)], Wmax)); R = cum(series([('f', 2, 2, n)], Wmax, pref=2**(n//2)))
        analyse(f"fermions n={n}: V=NS, M=NS", NS, NS, U, hl)
        analyse(f"fermions n={n}: V=NS, M=Ramond (ground 2^(n/2))", NS, R, U, hl)
        Rc = cum(series([('f', 2, 2, n)], Wmax, pref=2**(n//2-1)))
        analyse(f"fermions n={n}: V=NS, M=Ramond chiral ground 2^(n/2-1)", NS, Rc, U, hl)
    for n in (8, 32, 64):
        SB = cum(series([('b', 1, 2, 2*n)], Wmax)); analyse(f"symplectic bosons n={n} pairs (c=-n), M=V", SB, SB, U, hl)
        SF = cum(series([('f', 2, 2, 2*n)], Wmax)); analyse(f"symplectic fermions n={n} pairs (c=-2n), M=V", SF, SF, U, hl)
    # regraded complex fermions bc(D,1-D)^n, units 1/4
    U4 = 4; W4 = 2*Wmax
    for n in (8, 32, 64):
        NS = cum(series([('f', 2, 4, 2*n)], W4))        # 2n real NS fermions, weight 1/2 = 2 units
        for Dl in (1, 0):
            if Dl == 0: M = cum(series([('f', 4, 4, 2*n)], W4, pref=2**n))
            else: M = cum(series([('f', 1, 4, n), ('f', 3, 4, n)], W4))
            analyse(f"complex fermions n={n}: V=NS, M=bc({Dl}/4,1-{Dl}/4) regrading", NS, M, U4, list(range(2, W4//4, 2)))
