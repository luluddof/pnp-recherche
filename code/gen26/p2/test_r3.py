import numpy as np, math
rng = np.random.default_rng(5)
def mknet(g, kind):
    if kind == 'unif': c = rng.random((g, g))
    elif kind == 'hetero': u = rng.random(g)**4; c = np.outer(u, u)*rng.random((g, g))
    elif kind == 'star': c = np.zeros((g, g)); hub = rng.choice(g, max(3, g//8), replace=False); c[hub, :] = 1; c[:, hub] = 1
    elif kind == 'halfheavy': c = np.ones((g, g)); h = g//2; c[:h, :h] = 0; c[h:, h:] = 0   # bipartite
    elif kind == 'clique_half': c = np.zeros((g, g)); h = g//2; c[:h, :h] = 1
    c = np.triu(c, 1); c = c + c.T; np.fill_diagonal(c, 0); return c/c.max()
# (1) identity + drift exactly, with my own ANOVA
for kind in ('unif', 'hetero', 'star', 'halfheavy', 'clique_half'):
    g = 30; c = mknet(g, kind); U = c.sum(); mu = U/(g*(g-1)); ca = c.sum(1)
    al = (ca - (g-1)*mu)/(g-2); ga = c - mu - al[:, None] - al[None, :]; np.fill_diagonal(ga, 0)
    p = rng.permutation(g); e1 = 0; e2 = 0
    for j in range(0, g+1):
        S, T = p[:j], p[j:]
        F = c[np.ix_(S, T)].sum(); A = al[S].sum(); G = ga[np.ix_(S, T)].sum()
        e1 = max(e1, abs(F - (j*(g-j)*mu + (g-2*j)*A + G)))
        if j <= g-2:
            nxt = np.mean([ga[np.ix_(list(S)+[v], [w for w in T if w != v])].sum() for v in T])
            e2 = max(e2, abs(nxt - G*(1-2/(g-j))))
    # martingale increment bound |dX| <= 4 G_j / w_{j+1}
    w = lambda j: (g-j)*(g-j-1)/(g*(g-1)); worst = 0
    for j in range(0, g-2):
        S, T = list(p[:j]), list(p[j:]); G0 = ga[np.ix_(S, T)].sum(); Gm = max(abs(ga[v, S].sum()) for v in T) if j else 0
        for v in T:
            G1 = ga[np.ix_(S+[v], [x for x in T if x != v])].sum()
            dX = G1/w(j+1) - G0/w(j)
            if Gm > 1e-9: worst = max(worst, abs(dX)/(4*Gm/w(j+1)))
            elif abs(dX) > 1e-8: worst = float('inf')
    print(f"{kind:11s}: identity err {e1:.1e}, drift err {e2:.1e}, max |dX|/(4G/w) = {worst:.3f}")
# (2) proved bound (my own constants, R3 as stated, dA=dG=dGam=1/4 budget fix) vs random orders
def r3_bound_perN(g, nu, dA=0.25, dG=0.25, dGam=0.25):
    A = math.log(g/dA)*(g-1)/(nu*(g-2)**2)        # b7's (4x conservative) A-term
    R = 1/g + 2/(g-2); LG = math.log(g*g/dG); Lam = math.log(2*g/dGam)
    cum = 0; best = 0
    for j in range(1, g//2+1):
        i = j-1; wi1 = (g-i-1)*(g-i-2)/(g*(g-1)); cum += (4*R*math.sqrt(i*LG/2)/wi1)**2
        wj = (g-j)*(g-j-1)/(g*(g-1)); best = max(best, wj*math.sqrt(2*cum*Lam))
    return A, best/g
for g in (256, 1024, 4096):
    A, Gm = r3_bound_perN(g, 0.7071)
    print(f"g={g}: A-term/N {A:.4f}  Gamma-term/N {Gm:.4f}  (U/4 = nu/4 = 0.177)")
