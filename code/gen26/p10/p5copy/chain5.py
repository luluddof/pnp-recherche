# G26-P5: balanced-chain diagnostic of the overheads (greedy lt per level), model of P2/p9:
# e' = e (1/sqrt2 + z + 1/h) + 2^{-x}, x' = 2(x + lt), h = t K.  Decomposes ln(mu_L/mu_1) into growth / excess / 1/h / entries.
import math, sys, numpy as np
import mylib as M
Z = {'R2': M.zeta_R2, 'F': M.zeta_F, 'randfloor': M.zeta_randfloor, 'floor': M.zeta_floor, 'S2': M.zeta_S2, 'lin2': M.zeta_lin2}
LEV = 9
def chain(zn, lg, base='rowR3'):
    terms = (lambda l: (0.0, 0.0)) if base == 'rowzero' else M.r3_terms
    e = M.rowsweep_solver(lg, lg, terms); x = 2*lg
    mu1 = 1/(e*e*x); parts = np.zeros(4); rows = []
    for l in range(LEV):
        best = None
        for lt in np.arange(0.5, 48, 0.25):
            z = float(Z[zn](lt, x)); h = 2.0**(x + lt)
            e2 = e*(1/math.sqrt(2) + z + 1/h) + 2.0**-x; x2 = 2*(x + lt)
            mu2 = 1/(e2*e2*x2)
            if best is None or mu2 > best[0]: best = (mu2, lt, z, h, e2, x2)
        mu2, lt, z, h, e2, x2 = best
        parts += [math.log(x/(x+lt)), -2*math.log(1 + math.sqrt(2)*z), -2*math.log(1 + math.sqrt(2)/h/(1+math.sqrt(2)*z)),
                  -2*math.log(e2/(e*(1/math.sqrt(2) + z + 1/h)))]
        rows.append((x, lt, z*2.0**lt)); e, x = e2, x2
    mu = 1/(e*e*x)
    return mu1, mu, parts, rows, 2*lg
for zn in sys.argv[1:]:
    for base in ('rowR3', 'rowzero'):
        res = max((chain(zn, lg, base) for lg in np.arange(4, 30, 0.5)), key=lambda r: r[1])
        mu1, mu, parts, rows, x1 = res
        print(f"{zn:9s} {base:7s}: x1={x1:.0f} mu1={mu1:.3f} mu_L={mu:.4f} (c1 nosaw a=2: {2*math.sqrt(2/mu):.3f}); ln-loss growth {parts[0]:.3f} excess {parts[1]:.3f} 1/h {parts[2]:.4f} entries {parts[3]:.4f}")
        print("      levels (x, lt, z*t): " + " ".join(f"({r[0]:.0f},{r[1]:.2f},{r[2]:.1f})" for r in rows[:6]))
