# G26-P14: exact enumeration of stratified groupings (all (K!)^t), Var(Y12 - Y34) vs Var Y12 + Var Y34, and E|.| ratios.
import numpy as np, itertools, sys
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
def stats(c, K, t):
    P = np.array(list(itertools.permutations(range(K))))   # P[i, g] = node (within stratum) given to group g
    n = len(P)
    def Ypair(g, h):
        Y = np.zeros((n,)*t)
        for k in range(t):
            for l in range(t):
                C = c[k*K:(k+1)*K, l*K:(l+1)*K]
                if k == l:
                    v = C[P[:, g], P[:, h]]; shape = [1]*t; shape[k] = n; Y += v.reshape(shape)
                else:
                    M = C[P[:, g][:, None], P[:, h][None, :]]
                    shape = [1]*t; shape[k] = n; shape[l] = n
                    Y += M.reshape(shape) if k < l else M.T.reshape(shape) if False else np.moveaxis(M.reshape([n, n] + [1]*(t-2)), [0, 1], [k, l]) if t > 2 else (M if k < l else M.T)
        return Y.ravel()
    Y1 = Ypair(0, 1); Y2 = Ypair(2, 3)
    v1, v2 = Y1.var(), Y2.var(); vd = (Y1 - Y2).var()
    rv = vd/(v1 + v2) if v1 + v2 > 1e-12 else np.nan
    ra = np.abs(Y1 - Y2).mean()/(np.abs(Y1 - Y1.mean()).mean() + np.abs(Y2 - Y2.mean()).mean()) if v1 + v2 > 1e-12 else np.nan
    rs = np.abs(Y1 - Y2).mean()/(np.sqrt(2*max(v1, v2))) if v1 + v2 > 1e-12 else np.nan
    return rv, ra, rs, abs(Y1.mean() - Y2.mean())
def net(K, t, kind):
    g = K*t
    if kind == 'unif': c = rng.random((g, g))
    elif kind == 'bern': c = (rng.random((g, g)) < rng.random()).astype(float)
    elif kind == 'rank1': u = rng.random(g); c = np.outer(u, u)
    elif kind == 'star': c = np.zeros((g, g)); c[rng.integers(g), :] = 1; c[:, rng.integers(g)] = 1
    elif kind == 'block': lab = rng.integers(2, size=g); c = (lab[:, None] == lab[None, :]).astype(float)
    elif kind == 'sparse': c = (rng.random((g, g)) < 0.05).astype(float)
    c = np.triu(c, 1); c = c + c.T
    return c
if __name__ == '__main__':
    for K, t, reps in ((4, 2, 200), (4, 3, 100), (5, 2, 100), (5, 3, 20), (6, 2, 15)):
        eps = (4*K - 2)/((K - 2)*(K - 3))
        best = {}
        for kind in ('unif', 'bern', 'rank1', 'star', 'block', 'sparse'):
            mx = [0, 0, 0]; md = 0
            for r in range(reps):
                rv, ra, rs, dm = stats(net(K, t, kind), K, t)
                if np.isnan(rv): continue
                mx = [max(mx[0], rv), max(mx[1], ra), max(mx[2], rs)]; md = max(md, dm)
            best[kind] = mx
            print(f"K={K} t={t} {kind:6s}: max Var(Y-Y')/(VarY+VarY') = {mx[0]:.4f} (bound 1+eps={1+eps:.3f}); "
                  f"max E|Y-Y'|/(E|Y-m|+E|Y'-m|) = {mx[1]:.4f}; max E|Y-Y'|/(sqrt2 sigma) = {mx[2]:.4f} "
                  f"(lemma: {np.sqrt(1+eps):.3f}); |mean diff| max {md:.1e}", flush=True)
