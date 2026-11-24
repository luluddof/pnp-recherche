import glob, itertools
from core import *
# (a) q=1 construction W^T = W + sum_{a<b} tr(H_ab W) H_ab for N=2..7
for N in range(2, 8):
    col = [[None]*N for _ in range(N)]; sg = [[1]*N for _ in range(N)]; pal = {'I': 'S'}
    for a in range(N):
        for b in range(N):
            if a == b: col[a][b] = 'I'
            else: col[a][b] = ('H', min(a, b), max(a, b)); pal[col[a][b]] = 'T'; sg[a][b] = 1 if a < b else -1
    print('q=1 construction N', N, identity_tensor_ok(N, col, pal, sg), ranks(N, col, pal, sg), 'p expected', N*(N-1)//2)
# (b) balanced table
def pmin(N, q):
    p = 0
    while N*N > q*gstar(p, N) + p*gstar(q, N): p += 1
    return p
for N in (2, 3, 4, 5, 6, 7, 8, 16, 32, 64, 128):
    pb = next(p for p in range(1, N*N+1) if pmin(N, p) <= p)
    print('N', N, 'balanced counting bound', pb, ' N^(4/3)/2 =', round(N**(4/3)/2, 1))
# (c) HWY Prop 7.8 on the full submatrices of the U-arrays: each ordered colour pair in <= 2 positions
for f in sorted(glob.glob('/home/user/pnp-recherche/code/gen22/g22b/shape_*U16.json')):
    N, col, pal, sg = load_json(f)
    for k, cells in classes(N, col).items():
        s = len(cells); cnt = {}
        for i in range(s):
            for j1, j2 in itertools.combinations([j for j in range(s) if j > i], 2):
                key = (col[cells[i][0]][cells[j1][1]], col[cells[i][0]][cells[j2][1]])
                cnt[key] = cnt.get(key, 0) + 1
        mu = len(set(c for r in col for c in r))
        if s >= 3: print(f.split('/')[-1], 'U colour size', s, 'max positions per ordered pair', max(cnt.values()), 'C(s,3)', s*(s-1)*(s-2)//6, '<= 2 mu^2 =', 2*mu*mu)
