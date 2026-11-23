import json, sys, numpy as np
from collections import defaultdict
fn = sys.argv[1] if len(sys.argv) > 1 else 'results.jsonl'
nmin = int(sys.argv[2]) if len(sys.argv) > 2 else 50
V = ['base', 'bva', 'bva0', 'xor', 'and', 'orsub', 'rand']
D = defaultdict(lambda: defaultdict(list))
for line in open(fn):
    r = json.loads(line)
    for v in V:
        D[v][r['n']].append((r[v]['conf'], r[v]['proof'] + r[v]['defs']))
ns = sorted(D['base'])
print('median conflicts (median DRUP lines + def clauses)')
print('n    ' + ''.join(f'{v:>16}' for v in V))
for n in ns:
    print(f'{n:<5}' + ''.join(
        f"{int(np.median([c for c,_ in D[v][n]])):>8}/{int(np.median([p for _,p in D[v][n]])):<7}" for v in V),
        f' (k={len(D["base"][n])})')


def fit(x, y):
    A = np.vstack([x, np.ones_like(x)]).T
    co, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ co
    r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return co[0], r2

print(f'\nfits over n>={nmin}: exp: log2(med)=a*n+b -> base 2^a per var ; poly: log2(med)=k*log2(n)+b')
for metric, idx in (('conflicts', 0), ('proof', 1)):
    for v in V:
        xs = np.array([n for n in ns if n >= nmin], float)
        ys = np.array([np.log2(max(1, np.median([t[idx] for t in D[v][int(n)]]))) for n in xs])
        a, r2e = fit(xs, ys)
        k, r2p = fit(np.log2(xs), ys)
        # ratio to base
        print(f'{metric:9} {v:6} exp: 2^({a:.4f} n) = {2**a:.4f}^n R2={r2e:.4f} | poly: n^{k:.2f} R2={r2p:.4f} | better={"exp" if r2e>r2p else "poly"}')
# per-instance ratio variant/base
print('\nmedian per-instance ratio variant/base conflicts, by n')
for n in ns:
    print(n, ' '.join(f'{v}:{np.median([a[0]/max(1,b[0]) for a,b in zip(D[v][n],D["base"][n])]):.2f}' for v in V[1:]))
