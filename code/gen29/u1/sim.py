"""G29-U1 independent simulation of the unbalanced-level recursion (design level).
Coordinates: u = theta/2 (row cost per 2 lam), v = T/2. Frontier v = f(u), concave, f(0)=1, f(1)=0.
Level (my derivation from Prop 3 / Lemma 8(ii)):
  e' = max_u min( e_a (u + 1/h + 2 zeta_a) + 2^(1-x_a),  e_b (f(u) + 1/g + 2 zeta_b) + 2^(1-x_b) )
  x' = x_a + x_b + log2 t_a + log2 t_b,  t >= ln(2 g K)/(2 zeta^2)
Leaves: brute force e = 1/2 at every x >= 1.
DP over an integer grid x = 1..X (x = log2 K). Feasibility is monotone (t can be taken larger), so
e(x) = min(e(x-1), best split). Split candidates: x_b in a set S (all small values + geometric), x_a exact.
"""
import numpy as np, sys, math, time

def families():
    fam = {}
    fam['circle'] = lambda u: np.sqrt(np.clip(1 - u*u, 0, None))
    fam['line'] = lambda u: 1 - u
    fam['mono'] = lambda u: np.minimum(2*(1-u), 1 - u/2)          # p^2 q = N^4 (+ mirror)
    for c in (2, 4, 10):
        fam['circ_corner%d' % c] = (lambda c: lambda u: np.minimum.reduce(
            [np.sqrt(np.clip(1-u*u, 0, None)), c*(1-u), 1 - u/c]))(c)
    fam['circ_pow0.75'] = lambda u: np.minimum.reduce(
        [np.sqrt(np.clip(1-u*u, 0, None)), 2*(1-u)**0.75, 1 - (u/2)**(1/0.75)])
    return fam

def level(f, ea, eb, oa, ob, na, nb, it=60):
    """vectorised: max_u min(ea(u+oa)+na, eb(f(u)+ob)+nb); row increasing, col decreasing in u."""
    lo = np.zeros_like(ea); hi = np.ones_like(ea)
    for _ in range(it):
        m = 0.5*(lo+hi)
        r = ea*(m+oa)+na; c = eb*(f(m)+ob)+nb
        up = r < c
        lo = np.where(up, m, lo); hi = np.where(up, hi, m)
    m = 0.5*(lo+hi)
    return np.minimum(ea*(m+oa)+na, eb*(f(m)+ob)+nb)

def run(f, X, honest=True, Ssmall=200, ngeo=150):
    e = np.full(X+1, 0.5); e[0] = np.inf
    choice = np.zeros(X+1, dtype=np.int64)
    geo = np.unique(np.round(np.geomspace(1, X, ngeo)).astype(np.int64))
    S = np.unique(np.concatenate([np.arange(1, Ssmall+1), geo]))
    for x in range(2, X+1):
        if honest:
            Z = 2.0**(-np.arange(1, 31, 1.5))                    # zeta menu (same zeta both sides)
            T = (math.log(2)*(2*x+1))/(2*Z*Z)                      # t >= ln(2gK)/(2 zeta^2), ln(2gK) <= ln2(2x+1)
            D = 2*np.log2(T)
        else:
            Z = np.array([0.0]); T = np.array([1.0]); D = np.array([0.0])
        XB, ZZ = np.meshgrid(S, Z); TT = (math.log(2)*(2*x+1))/(2*ZZ*ZZ) if honest else np.ones_like(ZZ)
        DD = 2*np.log2(TT)
        XB = XB.ravel(); ZZ = ZZ.ravel(); TT = TT.ravel(); DD = DD.ravel()
        xa = np.floor(x - DD - XB).astype(np.int64)
        ok = (xa >= XB) & (XB >= 1)
        if not ok.any():
            e[x] = e[x-1]; continue
        xb = XB[ok]; xa = xa[ok]; zeta = ZZ[ok]; t = TT[ok]
        ea = e[xa]; eb = e[xb]
        if honest:
            oa = 2.0**(-xb.astype(float))/t + 2*zeta; ob = 2.0**(-xa.astype(float))/t + 2*zeta
            na = 2.0**(1-xa.astype(float)); nb = 2.0**(1-xb.astype(float))
        else:
            oa = ob = na = nb = np.zeros(len(xb))
        # both orientations (frontier may be asymmetric): big child on row side or on column side
        v1 = level(f, ea, eb, oa, ob, na, nb)
        v2 = level(f, eb, ea, ob, oa, nb, na)
        v = np.minimum(v1, v2)
        k = int(np.argmin(v))
        if v[k] < e[x-1]:
            e[x] = v[k]; choice[x] = xb[k]
        else:
            e[x] = e[x-1]; choice[x] = -1
    return e, choice

def local_alpha(e, xs):
    out = []
    for x1, x2 in zip(xs[:-1], xs[1:]):
        out.append((x1, x2, -math.log(e[x2]/e[x1])/math.log(x2/x1)))
    return out

if __name__ == '__main__':
    X = int(sys.argv[1]); names = sys.argv[2].split(',') ; honest = sys.argv[3] == 'honest'
    fam = families()
    for nm in names:
        t0 = time.time()
        e, ch = run(fam[nm], X, honest)
        xs = [x for x in (100, 300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000) if x <= X]
        la = local_alpha(e, xs)
        print('%s honest=%s X=%d (%.0fs)' % (nm, honest, X, time.time()-t0))
        for x in xs:
            print('   x=%7d e=%.4e x*e=%.3f sqrt(x)*e=%.3f choice x_b=%d' % (x, e[x], x*e[x], math.sqrt(x)*e[x], ch[x]))
        print('   local alpha:', ' '.join('[%d-%d] %.3f' % a for a in la))
        np.save('/tmp/claude-0/-home-user-pnp-recherche/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/gen29/u1/e_%s_%s_%d.npy' % (nm, 'h' if honest else 'c', X), e)
        sys.stdout.flush()
