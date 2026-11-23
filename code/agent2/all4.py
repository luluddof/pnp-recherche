"""Exact B2 circuit size of every 4-input Boolean function, via NPN classes (size is NPN-invariant in B2)."""
import numpy as np, itertools, json, time
from multiprocessing import Pool
from synth import circuit_size

n = 4; N = 16
bits = ((np.arange(65536)[:, None] >> np.arange(N)[None, :]) & 1).astype(np.uint8)
w = (1 << np.arange(N)).astype(np.int64)
canon = np.full(65536, 1 << 20, dtype=np.int64)
for perm in itertools.permutations(range(n)):
    for neg in range(N):
        # row t -> row t' where input perm[j] of t' = bit j of t xor neg
        src = []
        for t in range(N):
            tp = 0
            for j in range(n):
                tp |= (((t >> j) & 1) ^ ((neg >> j) & 1)) << perm[j]
            src.append(tp)
        g = bits[:, src] @ w
        canon = np.minimum(canon, g)
        canon = np.minimum(canon, (65535 - g))
reps, counts = np.unique(canon, return_counts=True)
print('NPN classes:', len(reps))

def job(f):
    s, exact = circuit_size(int(f), 4)
    return int(f), s, exact

if __name__ == '__main__':
    t0 = time.time()
    with Pool(4) as p:
        res = p.map(job, list(reps), chunksize=1)
    size = {f: s for f, s, e in res}
    assert all(e for _, _, e in res)
    dist = {}
    for f, c in zip(reps, counts):
        dist[size[int(f)]] = dist.get(size[int(f)], 0) + int(c)
    print('distribution over all 65536 functions:', dict(sorted(dist.items())))
    print('time', time.time() - t0)
    json.dump({'classes': {int(f): [size[int(f)], int(c)] for f, c in zip(reps, counts)}}, open('all4_sizes.json', 'w'))
