# Variants of chain_exact.py: (M) memoisation of tables (the additive 1/K_l entry term becomes a max, i.e. drops
# out of the main term), (O) optimistic row sweep with zs(g) = 1/(4(g-1)) + 2/g [CONJECTURAL ordering lemma].
import math, thmA_explicit as T, chain_exact as C
RHO = T.RHO
def chain(xs, es, z0, q, memo, L=26):
    out = [(xs, es)]; x, e = xs, es
    for l in range(L):
        z = z0*q**l; t = T.t_for(z, x); zs = T.zsum(t, x, True)
        g = t*2.0**x if x < 1000 else float('inf'); inv_g = 1/g if x < 1000 else 0.0
        invK = 0.0 if memo or x >= 1000 else 2.0**-x
        e_new = e*(RHO + inv_g + zs) + invK
        if memo: e_new = max(e_new, 2.0**-x*(1+RHO) if x < 1000 else 0)
        e = e_new; x = 2*(x + math.log2(t)); out.append((x, e))
    return out
for label, kind, memo, opt in [("clique base + memo", 'clique', True, False), ("row sweep + memo", 'row', True, False),
                               ("row sweep OPTIMISTIC ordering (conj.)", 'row', False, True)]:
    best = None
    starts = ([(float(x0), (((2**x0+1)**2)//4 - 1)/4**x0) for x0 in range(4, 16)] if kind == 'clique' else
              [(2*lg, (RHO + 2.0**-lg)/4 + ((1/(4*(2**lg-1)) + 2/2**lg) if opt else C.zs_sweep(lg))) for lg in [x/2 for x in range(12, 81)]])
    for xs, es in starts:
        for k in range(2, 12):
            for q in (0.4, 0.5, 0.6, 0.7, 0.8):
                ch = chain(xs, es, 2.0**-k, q, memo); mu = 1/(ch[-1][1]**2*ch[-1][0])
                if best is None or mu > best[0]: best = (mu, xs, es, k, q, ch)
    mu, xs, es, k, q, ch = best
    s2 = max(2*ch[l][1]*math.sqrt(2*ch[l+1][0]) for l in range(18, 25))
    print(f"{label:40s} mu_inf={mu:.4f} (x_start={xs}, e_start={es:.5f}, z0=2^-{k}, q={q}); single-schedule c1 a=2 {s2:.3f};"
          f" no-sawtooth a=2 {2*math.sqrt(2/mu):.3f} a=4 {2*math.sqrt(4/mu):.3f}")
