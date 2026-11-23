"""Exact check of the explicit degree-1 per-point construction (G6d), and comparison with the
common-E lower bound (G4c Hamming ball).
g_w(t) = chi(w) - t*chi(w)*r_w,  r_w = 1/tau_1 + (-1)^{w_1} * phi(w_{-1}),
phi(v) = (-1)^{|v|} * sum_{j>=2} v_j * (2/tau_j - 2/tau_1),   tau_k = k.
Claim: direction 1 perfect; every other direction has exactly 2^{l-2} bad edges => rank E_k = 2^{l-2}, rho=1/4.
For l<=7 also builds the 2^l x 2^l matrices J(tau_k) = H diag(g(tau_k)) H and checks
rank(J(tau_k) - K_k) with K_k = H (I_2 (x) D_+) H (exactly i-nice) equals the bad count."""
import itertools, sys
from fractions import Fraction as F
from math import comb
import sympy as sp

def construction(l):
    tau = [F(k) for k in range(1, l + 1)]
    V = list(itertools.product([0, 1], repeat=l))
    def r(w):
        v = w[1:]
        phi = (-1) ** sum(v) * sum(v[j - 1] * (2 / tau[j] - 2 / tau[0]) for j in range(1, l))
        return 1 / tau[0] + (-1) ** w[0] * phi
    g = {w: (lambda t, w=w: (-1) ** sum(w) * (1 - t * r(w))) for w in V}
    return tau, V, g

def bad_counts(l):
    tau, V, g = construction(l)
    per = [0] * l
    for w in V:
        for k in range(l):
            if w[k] == 0:
                u = list(w); u[k] = 1; u = tuple(u)
                if g[w](tau[k]) != g[u](tau[k]): per[k] += 1
    return per

def matrix_check(l):
    tau, V, g = construction(l)
    H1 = sp.Matrix([[1, 1], [1, -1]]); H = H1
    for _ in range(l - 1): H = sp.kronecker_product(H, H1)
    # index b in {0,1}^l (big-endian) <-> w ; diag entry g_w
    out = []
    for k in range(l):
        D = sp.diag(*[sp.Rational(g[w](tau[k]).numerator, g[w](tau[k]).denominator) for w in V])
        # nice part: replace g by its value with w_k = 0 (independent of w_k)
        Dn = sp.diag(*[sp.Rational(*(lambda x: (x.numerator, x.denominator))(g[tuple(w[:k]) + (0,) + tuple(w[k+1:])](tau[k]))) for w in V])
        E = H * (D - Dn) * H
        out.append(E.rank())
    return out

if __name__ == "__main__":
    for l in range(2, 11):
        per = bad_counts(l)
        lbcommon = sum(comb(l, j) for j in range(0, (l - 2) // 2 + 1))  # G4c common-E bound at d=1
        line = f"l={l}: bad per direction {per}  (2^(l-2)={2**(l-2)})  rho_max={max(per)/2**l}  common-E LB rho>={lbcommon/2**l:.3f}"
        if l <= 6: line += f"  matrix ranks {matrix_check(l)}"
        print(line, flush=True)
