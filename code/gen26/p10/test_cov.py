# exact check: (a,b) uniform on ordered distinct pairs of [K]; u,u' zero-mean vectors: E(u(a)+u'(b))^2 <= (1+1/(K-1))(su^2+su'^2)
# plus exact Doob-martingale check of R2F pieces on a tiny instance (K=3,t=3), all bijection tuples.
import numpy as np, itertools, math
rng = np.random.default_rng(1)
worst = 0
for K in (2, 3, 4, 5, 8):
    for _ in range(2000):
        u = rng.normal(size=K); u -= u.mean(); v = rng.normal(size=K) if rng.random() < .5 else -u*rng.random()+rng.normal(size=K)*.1; v -= v.mean()
        E = np.mean([(u[a]+v[b])**2 for a in range(K) for b in range(K) if a != b])
        bd = (1+1/(K-1))*(np.mean(u*u)+np.mean(v*v)); worst = max(worst, E/bd)
print("cov claim max ratio", worst)
# tiny exact R2F check
K, t = 3, 3; g = K*t; P = 1.0
for trial in range(5):
    c = rng.random((g, g))*P/g; c = np.triu(c, 1); c = c + c.T
    order = np.argsort(c.sum(1)); S = [order[k*K:(k+1)*K] for k in range(t)]
    perms = list(itertools.permutations(range(K)))
    # gamma=0, gamma'=1: a_k = S_k[pi^-1(0)], b_k = S_k[pi^-1(1)]; enumerate (a_k,b_k) ordered distinct pairs per stratum
    pairs = [[(S[k][i], S[k][j]) for i in range(K) for j in range(K) if i != j] for k in range(t)]
    def Y(sel): return sum(c[a, b] for (a, _) in sel for (_, b) in sel)
    allsel = list(itertools.product(*pairs))
    vals = {s: Y(s) for s in allsel}
    EY = np.mean(list(vals.values()))
    U = c.sum(); assert EY <= U/(K*(K-1)) + 1e-12
    # Doob martingale M_k
    def M(prefix):
        k = len(prefix); return np.mean([vals[s] for s in allsel if s[:k] == prefix])
    maxlin = 0; maxV = 0
    for s in allsel:
        lin = 0
        for k in range(t):
            Mk = M(s[:k+1]); Mk1 = M(s[:k])
            ca = c.sum(1); A = ca[S[k]].mean()
            lam = (ca[s[k][0]] - A)/K + (ca[s[k][1]] - A)/K
            lin += lam
        maxlin = max(maxlin, lin/(2*P/K))
    print("trial", trial, "sum Lambda / (2P/K) max over histories:", round(maxlin, 3))
