# Optimistic exponent for Casimir interpolation S = f(Omega_gl_n) on V = Sym^h(C^n) (bosons) or Lambda^h (fermions), n->inf.
# components of V(x)V: two-row (bosons) / two-column (fermions) shapes indexed by j; log-dim ~ n[g(eta+e)+g(eta-e)] for j=h-+e n
# interpolate the 2e n +1 components nearest the middle with degree d = 2 e n ; OSR >= dim(moves of <= d particles) >= n*2 g(d/n)
# (optimistic: ignores diagonal factors). p ~ remaining comps ~ n[g(eta+e)+g(eta-e)] ; N ~ n g(eta).
import numpy as np
H = lambda x: 0 if x <= 0 or x >= 1 else -(x*np.log2(x) + (1-x)*np.log2(1-x))
gB = lambda x: (1+x)*H(x/(1+x)) if x > 0 else 0.0
gF = lambda x: H(min(max(x, 0), 1))
for nm, g, etas in [('bosons', gB, np.linspace(0.02, 3, 150)), ('fermions', gF, np.linspace(0.02, 0.5, 49))]:
    best = (9, 0, 0)
    for eta in etas:
        for e in np.linspace(0, eta, 200):
            th = (g(eta+e) + g(eta-e)) / g(eta); T = 2*g(2*e) / g(eta)
            if nm == 'fermions' and eta + e > 1: continue
            m = max(th, T)
            if m < best[0]: best = (m, eta, e)
    print(nm, 'optimistic balanced exponent %.4f at eta=%.3f e=%.3f' % best)
best = (9, 0, 0)
for eta in np.linspace(0.02, 0.5, 49):
    for e in np.linspace(0, eta, 300):
        th = (gF(eta+e) + gF(eta-e)) / gF(eta); T = 2*max(gF(s) for s in np.linspace(0, e, 50)) / gF(eta)
        m = max(th, T)
        if m < best[0]: best = (m, eta, e)
print('fermions (fixed) optimistic balanced %.4f at eta=%.3f e=%.3f' % best)
