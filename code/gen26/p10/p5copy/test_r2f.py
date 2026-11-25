# G26-P5: checks for R2F. (1) exact conditional variance of Delta_k vs the bound (sqrt(2(1+1/(K-1))) (sd_Z + P/2g) + P/2g)^2
# along random histories; increments Delta_k <= 4 max|Z| + 3P/g; Lambda sum <= 2P/K. (2) stratified excess vs R2F bound.
import numpy as np, math, sys
sys.path.insert(0, '.')
import mylib as M
rng = np.random.default_rng(5)
def mknet(g, kind):
    if kind == 'unif': c = rng.random((g, g))
    elif kind == 'bern': c = (rng.random((g, g)) < 0.7)*1.0
    elif kind == 'hetero': u = rng.random(g)**4; c = np.outer(u, u)*rng.random((g, g))
    elif kind == 'rank2':  # c = (1 + e_a e_b)/2, e = +-1 balanced inside every sorted stratum (node sums equal)
        e = np.tile([1, -1], g//2); c = (1 + np.outer(e, e))/2
    elif kind == 'star': c = np.zeros((g, g)); hub = rng.choice(g, max(1, g//16), replace=False); c[hub, :] = 1; c[:, hub] = 1
    elif kind == 'blocks': lab = rng.integers(0, 2, g); c = (lab[:, None] == lab[None, :])*1.0
    c = np.triu(c, 1); c = c + c.T; np.fill_diagonal(c, 0)
    return c   # pair <= 1 =: P/g, P := g
def strata(c, K, t):
    order = np.argsort(c.sum(1), kind='stable'); return [order[i*K:(i+1)*K] for i in range(t)]
# (1) variance check
worst = 0; worstb = 0; worstl = 0
for kind in ('unif', 'hetero', 'rank2', 'star', 'blocks', 'bern'):
    K, t = 8, 12; g = K*t; c = mknet(g, kind); P = float(g); pg = 1.0
    S = strata(c, K, t); ca = c.sum(1)
    for rep in range(30):
        A = [rng.choice(S[k], 2, replace=False) for k in range(t)]   # (a_k, b_k) history
        lam = 0.0
        for k in range(t):
            Sk = S[k]; Sl = np.concatenate(S[k+1:]) if k+1 < t else np.array([], int)
            a_prev = np.array([A[i][0] for i in range(k)], int); b_prev = np.array([A[i][1] for i in range(k)], int)
            Z = np.array([sum(c[a, b_prev[j]] - c[a, S[j]].sum()/K for j in range(k)) for a in Sk])
            Zp = np.array([sum(c[a_prev[i], b] - c[S[i], b].sum()/K for i in range(k)) for b in Sk])
            f = Z - np.array([c[a, Sk].sum()/K for a in Sk]); fp = Zp - np.array([c[Sk, b].sum()/K for b in Sk])
            u = f - f.mean(); up = fp - fp.mean()
            Ew = c[np.ix_(Sk, Sk)].sum()/(K*(K-1))
            D = np.array([[u[i] + up[j] + c[Sk[i], Sk[j]] - Ew for j in range(K)] for i in range(K)])
            mask = ~np.eye(K, dtype=bool); vals = D[mask]
            assert abs(vals.mean()) < 1e-9
            var = (vals**2).mean()
            sdZ = max(Z.std(), Zp.std()); sd = max(np.sqrt((Z**2).mean()), np.sqrt((Zp**2).mean()))
            vb = (math.sqrt(2*(1+1/(K-1)))*(sd + pg/2) + pg/2)**2
            worst = max(worst, var/vb)
            zmax = max(np.abs(Z).max(), np.abs(Zp).max()); worstb = max(worstb, vals.max()/(4*zmax + 3*pg))
            # full D_k (b7 formula) minus Delta_k = Lambda_k
            lam += (ca[A[k][0]] - ca[Sk].mean())/K + (ca[A[k][1]] - ca[Sk].mean())/K
        worstl = max(worstl, lam/(2*P/K))
print(f"(1) max Var/bound {worst:.3f}  max Delta/(4 zmax+3P/g) {worstb:.3f}  max Lambda-sum/(2P/K) {worstl:.3f}  (all must be <= 1)")
# (2) stratified grouping: observed max excess vs R2F and R2 bounds (per N); rank2: compare with ln K / t
def Ymax(c, gr, K):
    G = np.zeros((len(gr), K)); G[np.arange(len(gr)), gr] = 1
    Y = G.T @ c @ G; np.fill_diagonal(Y, -np.inf); return Y.max()
for K, t in ((16, 64), (32, 32), (64, 16), (8, 256)):
    g = K*t; P = float(g); N = g*P; x = math.log2(K); lt = math.log2(t)
    bF = float(M.zeta_F(lt, x)); bR2 = float(M.zeta_R2(lt, x))
    for kind in ('unif', 'bern', 'hetero', 'rank2', 'star', 'blocks'):
        c = mknet(g, kind); U = c.sum(); S = strata(c, K, t); ex = []
        for r in range(8):
            gr = np.empty(g, int)
            for i in range(t): gr[S[i]] = rng.permutation(K)
            ex.append(K*K*Ymax(c, gr, K) - U)
        print(f"(2) K={K} t={t} {kind:6s}: max excess/N {max(ex)/N:+.4f} (x t = {max(ex)/N*t:+.2f})  R2F {bF:.4f}  R2 {bR2:.4f}  lnK/(2t)={math.log(K)/(2*t):.4f}", flush=True)
