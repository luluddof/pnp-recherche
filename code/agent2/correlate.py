"""Correlate candidate measures with exact B2 circuit size over all 65536 4-input functions (weighted by NPN class size)."""
import json, numpy as np, math
from measures import ALL
d = json.load(open('all4_sizes.json'))['classes']
reps = [int(k) for k in d]; size = np.array([d[str(f)][0] for f in reps]); w = np.array([d[str(f)][1] for f in reps])

def rankdata(x):
    o = np.argsort(x, kind='mergesort'); r = np.empty(len(x)); xs = x[o]; i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and xs[j + 1] == xs[i]: j += 1
        r[o[i:j + 1]] = (i + j) / 2; i = j + 1
    return r

def wcorr(a, b, w):
    a = a - np.average(a, weights=w); b = b - np.average(b, weights=w)
    return np.average(a * b, weights=w) / math.sqrt(np.average(a * a, weights=w) * np.average(b * b, weights=w))

# expand to all functions for Spearman (weights ok for Pearson; for Spearman expand)
print(f'{"measure":14s} pearson spearman  best-c s.t. c*mu<=C always   mean mu by size 0..7')
for name, m in ALL.items():
    mu = np.array([m([(f >> t) & 1 for t in range(16)]) for f in reps], dtype=float)
    pe = wcorr(mu, size.astype(float), w)
    ex_mu = np.repeat(mu, w); ex_s = np.repeat(size, w)
    sp = np.corrcoef(rankdata(ex_mu), rankdata(ex_s))[0, 1]
    ok = mu > 0
    c = (size[ok] / mu[ok]).min() if ok.any() else float('nan')
    means = [np.average(mu[size == s], weights=w[size == s]) for s in range(8)]
    print(f'{name:14s} {pe:6.3f}  {sp:6.3f}   c={c:5.3f}   ' + ' '.join(f'{x:5.2f}' for x in means))
