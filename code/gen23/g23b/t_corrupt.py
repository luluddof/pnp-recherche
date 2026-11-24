import glob, random, copy
from core import *
rnd = random.Random(7)
files = sorted(glob.glob('/home/user/pnp-recherche/code/gen22/g22b/shape_*.json')) + ['/home/user/pnp-recherche/code/gen22/g22b/grp_8_Q8.json'] + sorted(glob.glob('/home/user/pnp-recherche/code/gen22/g22b/transl_*.json'))
tot = mism = idtrue = 0
for f in files:
    N, col, pal, sg = load_json(f)
    for it in range(150):
        c2 = [r[:] for r in col]; s2 = [r[:] for r in sg]; p2 = dict(pal)
        kind = rnd.choice(['sign', 'pal', 'merge', 'swap', 'toU', 'recolour', 'allU', 'none'])
        a, b = rnd.randrange(N), rnd.randrange(N)
        ks = list(p2)
        if kind == 'sign': s2[a][b] *= -1
        elif kind == 'pal':
            k = rnd.choice(ks); p2[k] = rnd.choice([x for x in 'STU' if x != p2[k]])
        elif kind == 'merge':
            k1, k2 = rnd.sample(ks, 2); c2 = [[k1 if c == k2 else c for c in r] for r in c2]
        elif kind == 'swap':
            a2, b2 = rnd.randrange(N), rnd.randrange(N); c2[a][b], c2[a2][b2] = c2[a2][b2], c2[a][b]
        elif kind == 'toU':
            k = rnd.choice(ks); p2[k] = 'U'
        elif kind == 'recolour': c2[a][b] = ('new', 0); p2[('new', 0)] = rnd.choice('STU')
        elif kind == 'allU': p2 = {k: 'U' for k in p2}
        i1 = identity_tensor_ok(N, c2, p2, s2); pr = predicate(N, c2, p2, s2)
        tot += 1; idtrue += i1
        if i1 != pr: mism += 1; print('MISMATCH', f, kind)
print('corruption tests', tot, 'identity true', idtrue, 'mismatches', mism)
