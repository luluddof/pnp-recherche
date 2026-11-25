# G26-B7: explicit deviation terms of the new lemmas (see notes.md, sections R2, R3, R4).
# All quantities are per N (N = g P, P = max node sum, pair weights <= P/g), asymptotic model as in p9.
import math
import numpy as np
LN2 = math.log(2)
RHO = 1/math.sqrt(2)

# ---------- R2: stratified grouping (replaces Lemma 8 (ii) Hoeffding) ----------------------------------
# g = tK nodes -> K groups of exactly t nodes, strata of K consecutive nodes (sorted by node sum),
# independent uniform bijection per stratum.  For all gamma != gamma' (with prob >= 1 - d1 - d2):
#   K^2 Y <= U K/(K-1) + 4 N/t + N sqrt(2 Sig L2)/t^2,
#   Sig = sum_{k<t} (4 sqrt(k L/2) + 1)^2 <= 4 L t(t-1) + (16/3) sqrt(L/2) t^1.5 + t,
#   L = ln(2 g^2/d1),  L2 = ln(K(K-1)/(2 d2)).
def zeta_new(lt, x, nu=1.0, d1=0.5, d2=0.49):
    """lt = log2 t, x = log2 K (scalars or numpy arrays). Returns excess/N."""
    t = np.power(2.0, lt); lg = lt + x
    L = math.log(2/d1) + 2*lg*LN2
    L2 = 2*x*LN2 + math.log(1/(2*d2))           # ln(K^2/(2 d2)) >= ln(K(K-1)/(2 d2))
    Sig = 4*L*t*(t-1) + (16/3)*np.sqrt(L/2)*t**1.5 + t
    invK1 = np.where(x < 1000, 1/(np.power(2.0, np.minimum(x, 1000)) - 1), 0.0)
    return 4/t + np.sqrt(2*Sig*L2)/t**2 + nu*invK1

def zeta_old(lt, x):   # p2 split zeta (Hoeffding), = thmA_explicit.zsum(t, x, True)
    t = np.power(2.0, lt); lnK = x*LN2
    return np.sqrt((math.log(2.0001) + lnK)/(2*t)) + np.sqrt((np.log(2.0001*t) + 2*lnK)/(2*t))

# ---------- R3: ordering lemma for the row sweep (replaces zs0 of p9 Lemma RS) ---------------------------
# random order of g nodes; F_j = j(g-j) mu + (g-2j) A_j + Gamma_j (Hajek decomposition), padded to U = nu N.
#   max_j F_j <= U g/(4(g-1)) + max A^2/mu + max_j Gamma_j
#   A-term: ln(g/dA) (g-1)/(nu (g-2)^2) N
#   Gamma: X_j = Gamma_j/w_j martingale, |dX_i| <= 4 G_i / w_{i+1}, G_i = R sqrt(i LG/2), R = P/g + 2P/(g-2),
#          LG = ln(g^2/dG);  Gamma_j <= w_j sqrt(2 sum_{i<j} c_i^2 ln(2g/dGam)) for all j <= g/2 (and reverse).
def gamma_term(g, dG=1/3, dGam=1/3, npts=4000):
    """max_j Gamma_j / N for j <= g/2 (continuous approximation of the sums, exact w_j)."""
    R = 1/g + 2/(g-2)                         # in units of P
    LG = math.log(g*g/dG); Lam = math.log(2*g/dGam)
    # sum_{i<j} (4/w_{i+1})^2 R^2 i LG/2 ; w_i = (g-i)(g-i-1)/(g(g-1)) ~ (1-i/g)^2 (upper bound for 1/w)
    ys = np.linspace(0, 0.5, npts+1)
    f = 16*R*R*LG/2 * ys*g / ((1-ys-1/g)**2*(1-ys-2/g)**2) * (1.0)   # integrand in i, w_{i+1} lower bound
    cum = np.concatenate([[0], np.cumsum((f[1:]+f[:-1])/2*np.diff(ys)*g)])   # ~ sum over i < j
    w = (1-ys)*(1-ys-1/g)/(1-1/g)
    G = w*np.sqrt(2*cum*Lam)                  # in units of P
    return G.max()/g                          # per N (N = gP)

def zs_R3(lg, nu=RHO, dA=1/3):
    g = 2.0**lg
    A = math.log(g/dA)*(g-1)/(nu*(g-2)**2)
    return A + gamma_term(g) + 1/(4*(g-1)) + 2/g

def zs_p9(lg):
    g = 2.0**lg
    return math.sqrt(math.log(2*g*g + 2*g + 1)/(4*g))*(1 + 1/g) + 1/(4*(g-1)) + 2/g

if __name__ == "__main__":
    print("R2: zeta_old vs zeta_new at x=log K, t=2^lt")
    for x in (9, 20, 49, 137, 320):
        row = []
        for lt in (6, 8, 10, 12, 15, 19, 24):
            row.append(f"lt={lt}: {float(zeta_old(lt, x)):.2e}/{float(zeta_new(lt, x)):.2e}")
        print(f" x={x}: " + "  ".join(row))
    print("R3: zs_p9 vs zs_R3 (row sweep overhead per N)")
    for lg in (8, 10, 12, 14, 16, 16.5, 18, 20, 24):
        print(f" log g={lg}: p9 {zs_p9(lg):.3e}  R3 {zs_R3(lg):.3e}  (gamma {gamma_term(2.0**lg):.2e})")
