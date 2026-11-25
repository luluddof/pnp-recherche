import numpy as np, math
rng = np.random.default_rng(11)
def mknet(g, kind):
    if kind == 'unif': c = rng.random((g, g))
    elif kind == 'hetero': u = rng.random(g)**4; c = np.outer(u, u)*rng.random((g, g))
    elif kind == 'star': c = np.zeros((g, g)); hub = rng.choice(g, max(1, g//16), replace=False); c[hub, :] = 1; c[:, hub] = 1
    elif kind == 'ladder':  # node sums linearly spread, pairs concentrated near diagonal in sorted order
        idx = np.arange(g); c = (np.abs(idx[:, None]-idx[None, :]) < g//8)*(idx[:, None]+idx[None, :])/(2*g)
    elif kind == 'blocks':  # two communities, sorted interleaving
        lab = rng.integers(0, 2, g); c = (lab[:, None] == lab[None, :])*1.0
    c = np.triu(c, 1); c = c + c.T; np.fill_diagonal(c, 0)
    return c/c.max()   # pair weights <= 1 = P/g  => P = g
def Ymax(c, gr, K):
    g = len(gr); G = np.zeros((g, K)); G[np.arange(g), gr] = 1
    Y = G.T @ c @ G; np.fill_diagonal(Y, -np.inf); return Y.max()
for K, t in ((16, 64), (8, 256), (32, 32)):
    g = K*t; P = float(g); N = g*P
    L = math.log(2*g*g/0.5); L2 = math.log(K*(K-1)/(2*0.49))
    Sig = sum((4*math.sqrt(k*L/2)+1)**2 for k in range(t))
    bound = 4*N/t + N*math.sqrt(2*Sig*L2)/t**2
    for kind in ('unif', 'hetero', 'star', 'ladder', 'blocks'):
        c = mknet(g, kind); P_act = c.sum(1).max(); U = c.sum(); ca = c.sum(1); order = np.argsort(ca)
        # rescale so node sums <= P: here P_act may exceed g? pair <= 1, node sum <= g-1 < P. ok
        ex_s = []; ex_p = []
        for r in range(10):
            gr = np.empty(g, int)
            for i in range(t): gr[order[i*K:(i+1)*K]] = rng.permutation(K)
            ex_s.append(K*K*Ymax(c, gr, K) - U*K/(K-1))
            ex_p.append(K*K*Ymax(c, rng.permutation(np.repeat(np.arange(K), t)), K) - U*K/(K-1))
        print(f"K={K} t={t} {kind:7s}: max excess/N stratified {max(ex_s)/N:+.4f}  plain {max(ex_p)/N:+.4f}  R2 bound {bound/N:.4f} (lin part 4/t={4/t:.4f})")
