# Own implementation (G26-P2) of the deviation terms of R2, R3 and of the level formulas. Nothing imported from b7.
import math
import numpy as np
LN2 = math.log(2); RHO = 1/math.sqrt(2)

def zeta_R2(lt, x, d1=0.5, d2=0.49, myL=False):
    """R2 excess per N: 4/t + sqrt(2 Sig L2)/t^2 + 1/(K-1).  lt=log2 t, x=log2 K (numpy ok).
    L = ln(2 g^2/d1) as stated by b7 (myL: ln(4 g K/d1), my union count 2gK events x 2 tails)."""
    lt = np.asarray(lt, float); x = np.asarray(x, float)
    t = np.power(2.0, lt)
    L = (math.log(4/d1) + (2*x + lt)*LN2) if myL else (math.log(2/d1) + 2*(lt + x)*LN2)
    # L2 = ln(K(K-1)/(2 d2)) <= 2 x ln2 - ln(2 d2)
    L2 = 2*x*LN2 - math.log(2*d2)
    Sig = 4*L*t*(t-1) + 8*np.sqrt(L/2)*(2/3)*t**1.5 + t
    invK1 = np.where(x < 60, 1/(np.power(2.0, np.minimum(x, 60)) - 1), np.power(2.0, -np.minimum(x, 1e4)))
    return 4/t + np.sqrt(2*Sig*L2)/t**2 + invK1

def zeta_hoeff(lt, x):   # p2 split-zeta Hoeffding (old), own transcription of the formula in p2 notes
    t = np.power(2.0, np.asarray(lt, float)); lnK = np.asarray(x, float)*LN2
    return np.sqrt((math.log(2.0001) + lnK)/(2*t)) + np.sqrt((np.log(2.0001*t) + 2*lnK)/(2*t))

def r3_terms(lg, dA=0.33, dG=0.33, dGam=0.33, b7const=False, npts=20000):
    """R3 (random order) deviation terms per N for a network on g = 2^lg nodes, P = N/g, pair <= P/g.
    returns (cA, Gam): max_j F_j <= U g/(4(g-1)) + cA*N^2/U_pad + Gam*N  (cA: A^2/mu = cA N / nu).
    Own constants: A-range P/(g-2); gamma-range P/g + P/(g-2); G-events forward+reverse (g^2, two-sided);
    Azuma one-sided, g events (j <= g/2, two directions).  b7const=True: b7's constants and budget."""
    g = 2.0**lg
    if b7const:
        cA = math.log(g/(1/3))*(g-1)/((g-2)**2)
        R = 1/g + 2/(g-2); LG = math.log(g*g/(1/3)); Lam = math.log(2*g/(1/3))
    else:
        cA = math.log(g/dA)*(g-1)/(4*(g-2)**2)
        R = 1/g + 1/(g-2); LG = math.log(2*g*g/dG); Lam = math.log(g/dGam)
    # S(j) = sum_{i<j} c_i^2, c_i = 4 R sqrt(i LG/2)/w_{i+1} (units of P), 1/w_{i+1} <= 1/(1-(i+2)/g)^2
    # <= integral_0^j of the (increasing, convex) bound; trapezoid overestimates convex integrals.
    ys = np.linspace(0, 0.5, npts+1); i = ys*g
    f = 16*R*R*(LG/2)*i/np.maximum(1-(i+2)/g, 1e-300)**4
    S = np.concatenate([[0], np.cumsum((f[1:]+f[:-1])/2*np.diff(i))])
    w = (g - i)*(g - i - 1)/(g*(g-1))
    Gam = (w*np.sqrt(2*S*Lam)).max()/g
    return cA, Gam

def rowsweep_cost(xr, lg, cA, Gam):
    """per-N cost of sweeping a network with U/N = xr on g=2^lg nodes (order from R3), padding optimised."""
    g = 2.0**lg
    nus = 2*math.sqrt(cA)
    v = np.maximum(xr, nus)
    return xr*g/(4*(g-1)) + (v - xr)/4 + cA/v + Gam + 2/g

def rowsweep_solver(lg, lh, terms):
    """g x h grid (x' = lg + lh), sweep rows (g nodes, U_row/N = 1/h + cos th) or columns (h nodes, 1/g + sin th),
    whichever is cheaper per iota; worst iota = max over the arc (radius 1). Bisection on the crossing."""
    cAg, Gg = terms(lg); cAh, Gh = terms(lh)
    g, h = 2.0**lg, 2.0**lh
    X = lambda th: rowsweep_cost(1/h + math.cos(th), lg, cAg, Gg)   # nonincreasing in th
    Y = lambda th: rowsweep_cost(1/g + math.sin(th), lh, cAh, Gh)   # nondecreasing
    if X(0) <= Y(0): return X(0)
    if X(math.pi/2) >= Y(math.pi/2): return Y(math.pi/2)
    lo, hi = 0.0, math.pi/2
    for _ in range(60):
        mid = (lo+hi)/2
        if X(mid) > Y(mid): lo = mid
        else: hi = mid
    return max(min(X(lo), Y(lo)), min(X(hi), Y(hi)))

def level_U(ea, eb, al, be):
    """max over the arc of min(ea (cos th + al), eb (sin th + be)), vectorised (closed form)."""
    ea, eb, al, be = np.broadcast_arrays(np.asarray(ea, float), np.asarray(eb, float), np.asarray(al, float), np.asarray(be, float))
    X0 = ea*(1+al); Y0 = eb*be; X1 = ea*al; Y1 = eb*(1+be)
    R = np.sqrt(ea*ea + eb*eb); phi = np.arctan2(eb, ea)
    c = np.clip((eb*be - ea*al)/R, -1, 1)
    th = np.arccos(c) - phi
    th = np.clip(th, 0, math.pi/2)
    val = np.minimum(ea*(np.cos(th)+al), eb*(np.sin(th)+be))
    val = np.where(X0 <= Y0, X0, np.where(X1 >= Y1, Y1, val))
    return val

def s_sweep(K): return ((K+1)**2)//4 - 1
def d_dc(K):   # balanced D&C on the middle cut: d(K) = floor(K/2) ceil(K/2) + d(ceil(K/2))
    return 0 if K <= 1 else (K//2)*(K - K//2) + d_dc(K - K//2)

# ---------------- G26-P5 additions ----------------
def zeta_F(lt, x, d1=0.33, d2=0.33, d3=0.33, mcd=True):
    """R2F (Freedman version of R2), excess per N:  2/t + 1/(K-1) + s'/t^2,
    s' = b'L2/3 + sqrt((b'L2/3)^2 + 2 V' L2)   (Freedman, one-sided, union over K(K-1)/2 pairs),
    V' = sum_{k=1..t} (s2 (m sqrt(k-1) + 1/2) + 1/2)^2,  s2 = sqrt(2(1+1/(K-1))),
    m = min(sqrt(L/2), 1/2 + sqrt(ell/2))  [Hoeffding sup bound / McDiarmid rms bound],
    b' = 4 sqrt((t-1)L/2) + 3,  L = ln(4gK/d1), ell = ln(2Kt/d3), L2 = ln(K(K-1)/(2 d2))."""
    lt = np.asarray(lt, float); x = np.asarray(x, float)
    t = np.power(2.0, lt)
    L = math.log(4/d1) + (2*x + lt)*LN2
    ell = math.log(2/d3) + (x + lt)*LN2
    L2 = 2*x*LN2 - math.log(2*d2)
    Km1 = np.where(x < 60, np.power(2.0, np.minimum(x, 60)) - 1, np.inf)
    s2 = np.sqrt(2*(1 + 1/Km1))
    m = np.minimum(np.sqrt(L/2), 0.5 + np.sqrt(ell/2)) if mcd else np.sqrt(L/2)
    V = s2*s2*m*m*t*(t-1)/2 + 2*s2*m*(s2+1)/2*(2/3)*t**1.5 + t*((s2+1)/2)**2
    b = 4*np.sqrt(np.maximum(t-1, 0)*L/2) + 3
    s = b*L2/3 + np.sqrt((b*L2/3)**2 + 2*V*L2)
    invK1 = np.where(x < 60, 1/(np.power(2.0, np.minimum(x, 60)) - 1), np.power(2.0, -np.minimum(x, 1e4)))
    return 2/t + s/t**2 + invK1

from scipy.special import ndtri, log_ndtr
from scipy.optimize import brentq
def zeta_floor(lt, x, nu=1/math.sqrt(2)):
    """Counting (first-moment) floor: random network c_ab = (P/g) Bern(nu); every grouping (any group sizes; at most K^g
    groupings) has a pair with K^2 Y >= U + kappa y* N / t,  kappa = sqrt(nu(1-nu)),  Phibar(y*) = 2 g ln K/(K(K-1))
    (normal approximation of the binomial tail)."""
    lt = np.asarray(lt, float); x = np.asarray(x, float)
    lt, x = np.broadcast_arrays(lt, x); sh = lt.shape
    out = np.zeros(lt.size)
    for i, (a, b) in enumerate(zip(lt.ravel(), x.ravel())):
        lnK = b*LN2
        lq = math.log(2*lnK) + a*LN2 - (math.log(2.0**b - 1) if b < 60 else lnK)
        if lq >= math.log(0.5): y = 0.0
        else: y = brentq(lambda yy: log_ndtr(-yy) - lq, 0.0, 2000.0)
        out[i] = math.sqrt(nu*(1-nu))*y/2.0**a
    return out.reshape(sh) if sh else float(out[0])

def zeta_lin2(lt, x):   # linear part only of R2F (degenerate part free): 2/t + 1/(K-1)
    t = np.power(2.0, np.asarray(lt, float)); x = np.asarray(x, float)
    return 2/t + np.where(x < 60, 1/(np.power(2.0, np.minimum(x, 60)) - 1), 0.0)

def zeta_randfloor(lt, x):   # behaviour of random stratified grouping on the rank-2 network: excess ~ ln K / t (+0 linear)
    return np.asarray(x, float)*LN2/np.power(2.0, np.asarray(lt, float))

def zeta_S2(lt, x):
    """Structural option S2 (edge-weighted cliques, per-rectangle padding): effective excess per N
    = (1/2) K^2 max_pair sd(Y)/N + 1/(K-1),  K^2 sd(Y)/N <= (1/t)[ 1/sqrt2 + sqrt(sum_k (s2 (sqrt(k-1)/2 + 1/2) + 1/2)^2)/t ]."""
    lt = np.asarray(lt, float); x = np.asarray(x, float); t = np.power(2.0, lt)
    Km1 = np.where(x < 60, np.power(2.0, np.minimum(x, 60)) - 1, np.inf)
    s2 = np.sqrt(2*(1 + 1/Km1))
    V = s2*s2*t*(t-1)/8 + 2*s2*0.5*(s2+1)/2*(2/3)*t**1.5 + t*((s2+1)/2)**2
    sd = (1/math.sqrt(2) + np.sqrt(V)/t)/t
    return 0.5*sd + np.where(x < 60, 1/Km1, 0.0)

def zeta_zero(lt, x):   # no grouping excess at all (only the K/(K-1) term): shows the remaining non-grouping overheads
    x = np.asarray(x, float) + 0*np.asarray(lt, float)
    return np.where(x < 60, 1/(np.power(2.0, np.minimum(x, 60)) - 1), 0.0)

def r3_terms_v2(lg, dA=0.33, dG=0.33, dGam=0.33, npts=20000):
    """R3'' (G26-P5): P2's R3 terms with (i) Hoeffding-Azuma in its interval form for the Gamma martingale (the increment
    dX_j = -2(gamma(v,S_j) - avg)/w_{j+1}, v uniform on S_j^c, lies in a predictable interval of length <= 4G_j/w_{j+1}):
    deviation sqrt(S Lam/2) instead of sqrt(2 S Lam); (ii) Serfling's finite-population factor (1 - (i-1)/(g-1)) in the
    sup bound G_i (gamma(v, S_i) = sum over a uniform i-subset of [g]-{v}) and (1 - (j-1)/g) in the A-term (max at j = g/2)."""
    g = 2.0**lg
    cA = math.log(g/dA)*(g-1)/(4*(g-2)**2)*(0.5 + 1/g)
    R = 1/g + 1/(g-2); LG = math.log(2*g*g/dG); Lam = math.log(g/dGam)
    ys = np.linspace(0, 0.5, npts+1); i = ys*g
    f = 16*R*R*(LG/2)*i*np.maximum(1 - (i-1)/(g-1), 0)/np.maximum(1-(i+2)/g, 1e-300)**4
    # f is no longer convex near 0 only through the (i-1) shift; trapezoid on 20000 points: error negligible, add 0.1% margin
    S = np.concatenate([[0], np.cumsum((f[1:]+f[:-1])/2*np.diff(i))])*1.001
    w = (g - i)*(g - i - 1)/(g*(g-1))
    Gam = (w*np.sqrt(S*Lam/2)).max()/g
    return cA, Gam
