import json, numpy as np, sys
def load(t):
    d = json.load(open(f'sol_{t}.json')); X = np.array(d['x']); E = np.array(d['e']); o = np.argsort(X)
    return X[o], np.minimum.accumulate(E[o])
Xb, Eb = load('cinf')
grid = np.exp(np.linspace(np.log(8), np.log(30000), 4000))
def ev(X, E, g): k = np.searchsorted(X, g, side='right') - 1; return E[np.maximum(k, 0)]
eb = ev(Xb, Eb, grid)
for t in sys.argv[1:]:
    X, E = load(t); ec = ev(X, E, grid); r = eb/ec
    out = []
    for thr in (1.01, 1.1, 1.5, 2.0, 4.0):
        i = np.argmax(r >= thr)
        out.append(f"gain>={thr}: x={grid[i]:.0f}" if r[i] >= thr else f"gain>={thr}: none")
    pts = ' '.join(f"{x}:{eb[np.searchsorted(grid,x)]/ec[np.searchsorted(grid,x)]:.2f}" for x in (32,64,128,256,512,1000,2000,4000,8000,16000))
    print(t, '|', '; '.join(out), '| ratio e_bal/e_c at x:', pts)
