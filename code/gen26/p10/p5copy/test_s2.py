# G26-P5: Lemma S2 check: empirical sd of Y_{gg'} (stratified random grouping) and mean |Y - EY| vs the bound
# sd <= (P/K)[1/sqrt2 + sqrt(sum_k (s2(sqrt(k-1)/2+1/2)+1/2)^2)/t]; also padding cost of a random placement vs (1/2) sum|Y-EY|.
import numpy as np, math
rng = np.random.default_rng(7)
def mknet(g, kind):
    if kind == 'bern': c = (rng.random((g, g)) < 0.7)*1.0
    elif kind == 'unif': c = rng.random((g, g))
    elif kind == 'rank2': e = np.tile([1, -1], g//2); c = (1 + np.outer(e, e))/2
    elif kind == 'jump': c = np.ones((g, g)); light = rng.choice(g, g//(2*16) if g >= 32 else 1, replace=False); c[light, :] = 0; c[:, light] = 0
    elif kind == 'hetero': u = rng.random(g)**4; c = np.outer(u, u)*rng.random((g, g))
    elif kind == 'star': c = np.zeros((g, g)); hub = rng.choice(g, max(1, g//16), replace=False); c[hub, :] = 1; c[:, hub] = 1
    c = np.triu(c, 1); c = c + c.T; np.fill_diagonal(c, 0); return c
for K, t in ((16, 16), (16, 64), (32, 32)):
    g = K*t; P = float(g)
    s2 = math.sqrt(2*(1+1/(K-1))); V = sum((s2*(math.sqrt(k)/2 + 0.5) + 0.5)**2 for k in range(t))
    bound = (P/K)*(1/math.sqrt(2) + math.sqrt(V)/t)
    for kind in ('bern', 'unif', 'rank2', 'jump', 'hetero', 'star'):
        c = mknet(g, kind); order = np.argsort(c.sum(1), kind='stable'); S = [order[i*K:(i+1)*K] for i in range(t)]
        Ys = []
        for r in range(60):
            gr = np.empty(g, int)
            for i in range(t): gr[S[i]] = rng.permutation(K)
            G = np.zeros((g, K)); G[np.arange(g), gr] = 1; Y = G.T @ c @ G
            Ys.append(Y[np.triu_indices(K, 1)])
        Ys = np.array(Ys); EY = Ys.mean()
        sd = Ys.std(); mad = np.abs(Ys - EY).mean()
        print(f"K={K} t={t} {kind:6s}: sd/bound {sd/bound:.3f}  mean|Y-EY|/bound {mad/bound:.3f}  K^2 mad/(2N) * t = {K*K*mad/(2*g*P)*t:.3f} (S2 constant 0.60)")
