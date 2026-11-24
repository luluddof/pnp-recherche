# Optimal valuation taper for one App.-D level, BOTH counting regimes (p>=n partition regime, p<n binomial regime).
# Units: M = 1.  Levels i=0..u-1, D_i = 2^(a-i), n_i = nu/D_i modes per weight, p_i = D_i*Lambda_i (+ evalloss).
# log #tuples ~ f(p,n) = min_theta theta*p + n*sum_j -ln(1-exp(-theta j))   (saddle-point; exact to leading order).
import numpy as np, math, sys
from scipy.optimize import linprog, minimize_scalar, brentq

def g(theta):
    if theta < 0.05:   # Dedekind-eta asymptotics, error exp(-4pi^2/theta)
        return math.pi**2/(6*theta) + 0.5*math.log(theta/(2*math.pi)) - theta/24
    j = np.arange(1, int(60/theta)+2)
    return float(np.sum(-np.log1p(-np.exp(-theta*j))))
def f(p, n):
    if p <= 0: return 0.0
    res = minimize_scalar(lambda lt: math.exp(lt)*p + n*g(math.exp(lt)), bounds=(-40, 6), method='bounded',
                          options={'xatol':1e-10})
    return res.fun
def lam_max(E, D, nu):
    n = nu/D
    h = lambda lam: f(D*lam, n) - E
    lo, hi = 0.0, 1.0
    while h(hi) < 0:
        hi *= 2
        if hi > 1e9: return hi
    return brentq(h, lo, hi, xtol=1e-12)
def feasible(E, a, u, nu):
    lm = np.array([lam_max(E, 2.0**(a-i), nu) for i in range(u)])
    A = np.zeros((u, u))
    for i in range(u):
        for j in range(1, a-i+1):
            w = j/2.0**j
            A[i, i] += w
            if i+j < u: A[i, i+j] -= w
    c = -np.array([2.0**(-i) for i in range(u)])
    r = linprog(c, A_ub=A, b_ub=lm, bounds=[(0, None)]*u, method='highs')
    return (-r.fun) >= 2 - 2.0**(-a), r.x, lm
def best(a, u, nu):
    lo, hi = 0.0, 50.0
    for _ in range(40):
        mid = (lo+hi)/2
        ok, q, lm = feasible(mid, a, u, nu)
        if ok: hi = mid
        else: lo = mid
    ok, q, lm = feasible(hi, a, u, nu)
    return hi, q, lm
print("# E* = min over tapers of max_i ln(#tuples_i)/M ;  report E*, E*sqrt(u), and paper-style bound")
for nu in (2.0, 6.0):
    for (a, u) in [(4,4),(8,8),(12,12),(16,16),(24,24),(32,32),(48,48),(16,8),(32,16),(48,24),(60,30),(20,16),(40,32)]:
        E, q, lm = best(a, u, nu)
        print(f"nu={nu} a={a:3d} u={u:3d}  E*={E:.4f}  E*sqrt(u)={E*math.sqrt(u):.4f}  E*sqrt(u/nu)={E*math.sqrt(u/nu):.4f}  q0..={np.round(q[:4],3)} q_last={np.round(q[-3:],4)}")
        sys.stdout.flush()
