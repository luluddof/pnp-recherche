"""Exp A2: rank of the product C of crossing clauses alone, as a function of k = |dA| = |dB| and crossing density g = c/k.
Each crossing clause: 3 distinct variables, split 1+2 or 2+1 between the sides at random, random signs.
Every boundary variable is forced to occur (otherwise it would not be a boundary variable): we resample until all occur.
Reports log2(rank C)/k (rank over R via SVD; 0/1-ish integer matrices, sizes <= 512)."""
import numpy as np, random, sys
rng = random.Random(5)
def inst(k, c):
    while True:
        cls = []
        for _ in range(c):
            na = rng.choice([1, 2])
            a = rng.sample(range(k), na); b = rng.sample(range(k), 3 - na)
            cls.append(([(v, rng.random() < .5) for v in a], [(v, rng.random() < .5) for v in b]))
        if len({v for ca, _ in cls for v, _ in ca}) == k and len({v for _, cb in cls for v, _ in cb}) == k: return cls
X = lambda k: (np.arange(2 ** k)[:, None] >> np.arange(k)[None, :]) & 1
print('k g  mean log2rank/k  [min,max] over 4 samples;  mean log2(#distinct row patterns)/k')
for g in [1.0, 2.0, 4.0]:
    for k in [4, 6, 8, 9]:
        c = max(k, int(g * k)); xa = X(k); res = []; pats = []
        for rep in range(4):
            cls = inst(k, c)
            C = np.ones((2 ** k, 2 ** k))
            F = []
            for ca, cb in cls:
                fa = np.all([xa[:, v] == neg for v, neg in ca], axis=0); fb = np.all([xa[:, v] == neg for v, neg in cb], axis=0)
                C *= 1 - np.outer(fa, fb); F.append(fa)
            res.append(np.log2(max(1, np.linalg.matrix_rank(C))) / k)
            pats.append(np.log2(len({tuple(r) for r in np.array(F).T})) / k)
        print(k, g, f'{np.mean(res):.3f} [{min(res):.3f},{max(res):.3f}]  {np.mean(pats):.3f}', flush=True)
