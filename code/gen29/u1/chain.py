"""Continuous big-side chain (ODE in x) for corner families, small children from the DP table (honest overheads).
dPsi/dx = max_{x_b, zeta} (1/e_new - 1/e)/(x_b + delta)  where e_new = level(e, e_table[x_b]) at x + x_b + delta.
Achievable (lower bound on Psi): each step is a legitimate Lemma-U level with the current solver as big child."""
import numpy as np, math, sys
sys.path.insert(0, '/tmp/claude-0/-home-user-pnp-recherche/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/gen29/u1')
from sim import families, level
D = '/tmp/claude-0/-home-user-pnp-recherche/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/gen29/u1/'
nm = sys.argv[1]; X1 = int(sys.argv[2]); xend = float(sys.argv[3]); bf = len(sys.argv) > 4 and sys.argv[4] == 'bf'
f = families()[nm]
tab = np.load(D + 'e_%s_h_%d.npy' % (nm, X1))
if bf: tab = np.full_like(tab, 0.5)            # small children = brute force only
S = np.unique(np.round(np.geomspace(1, X1, 400)).astype(np.int64))
Zm = 2.0**(-np.arange(1, 60, 1.0))
XB, ZZ = np.meshgrid(S, Zm); XB = XB.ravel(); ZZ = ZZ.ravel()
eb = tab[XB]
x = float(X1); e = tab[X1]
eta = 2e-3; out = []
nxt = 10.0**math.ceil(math.log10(x))
while x < xend:
    TT = math.log(2)*(2*x+1)/(2*ZZ*ZZ); dl = 2*np.log2(TT)
    oa = 2.0**(-XB.astype(float))/TT + 2*ZZ; ob = 2.0**(-x)/TT + 2*ZZ
    na = 2.0**(1-x); nb = 2.0**(1-XB.astype(float))
    ea = np.full(len(XB), e)
    en = np.minimum(level(f, ea, eb, oa, ob, na, nb), level(f, eb, ea, ob, oa, nb, na))
    r = (1/en - 1/e)/(XB + dl)
    r = np.where(XB + dl <= 0.05*x, r, -np.inf)
    k = int(np.argmax(r))
    dx = eta*x
    Psi = 1/e + r[k]*dx; x += dx; e = 1/Psi
    if x >= nxt:
        out.append((x, e, XB[k], ZZ[k], dl[k]))
        print('x=%.3e  x*e=%.3f  x*e/log2x=%.4f  x*e/(log2 x)^0.5=%.3f  chosen x_b=%d zeta=%.1e delta=%.1f' % (
            x, x*e, x*e/math.log2(x), x*e/math.sqrt(math.log2(x)), XB[k], ZZ[k], dl[k]))
        sys.stdout.flush()
        nxt *= 10
