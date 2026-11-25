import json, numpy as np
from spread import swapT
B = '/home/user/pnp-recherche/code/'
def load_g19(f):
    d = json.load(open(B + 'gen19/g19b/' + f)); N = int(d['N'])
    G = np.array(d['G'], dtype=float); Gp = np.array(d['Gp'], dtype=float)
    C = np.einsum('lxu,lyv->xyuv', G, Gp); return swapT(N) - C, C, int(d['p']), int(d['q'])
def load_g24g():
    d = json.load(open(B + 'gen24/g24g/exact_6_12_3.json')); N = int(d['N'])
    E = np.array(d['E'], dtype=float); Cc = np.array(d['C'], dtype=float)
    C = np.einsum('lxu,lyv->xyuv', E, Cc); return swapT(N) - C, C, int(d['rank']), int(d['osr'])
def load_mono(f):   # g22b STU arrays (U empty), as in gen25/b2/tdual.py
    d = json.load(open(f)); N = d['N']; cells = d['cells']; cls = {}
    for a in range(N):
        for b in range(N):
            pal, k, s = cells[a][b]; cls.setdefault((pal, k), np.zeros((N, N)))[a, b] = s
    C = np.zeros((N,)*4); R = np.zeros((N,)*4)
    for (pal, k), K in cls.items():
        if pal == 'S': C += np.einsum('xu,vy->xyuv', K, K)
        else: R += np.einsum('vu,xy->xyuv', K, K)
    return R, C, d['p'], d['q']
def appb(h, N=None):
    import appb_exact as ab, math
    N, mus, I, T = ab.pieces(h, N); d = math.factorial(h)
    A = {s: np.array(T[s], dtype=float) / d for s in I}
    return N, I, A
def hrubes(r, k):
    import hr
    S = hr.setup(r, k); tm = hr.terms(S); N = S['N']
    C = sum(np.einsum('xu,yv->xyuv', G.astype(float), Gp.astype(float)) for G, Gp in tm) / (2 * S['den'])
    return swapT(N) - C, C, S
def cd(n):
    from cd import table
    S, K = table(n); N = 1 << n; a = np.arange(N); C = np.zeros((N,)*4)
    for x in range(N):
        c = a ^ x; G = np.zeros((N, N)); Gp = np.zeros((N, N))
        G[a, c] = S[a, c]; Gp[a, c] = S[c, a]
        C += np.einsum('xu,yv->xyuv', G, Gp) / 2
    return swapT(N) - C, C
