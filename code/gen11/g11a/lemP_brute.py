# brute force: for random 3-CNF, n small, Delta=2: is there S with >= m/6 inside, >= m/6 inside complement, <= m/20 crossing?
import numpy as np, random, sys
def run(n, Delta, seed):
    rnd = random.Random(seed); m = int(Delta*n)
    masks = np.array([sum(1 << v for v in rnd.sample(range(n), 3)) for _ in range(m)], dtype=np.int64)
    S = np.arange(1 << (n-1), dtype=np.int64)   # S and complement symmetric: fix vertex n-1 outside S
    full = (1 << n) - 1
    inside = np.zeros(S.shape, dtype=np.int16); outside = np.zeros(S.shape, dtype=np.int16)
    for mk in masks:
        inside += ((S & mk) == mk)
        outside += ((S & mk) == 0)
    cross = m - inside - outside
    bad = (6*inside >= m) & (6*outside >= m) & (20*cross <= m)
    # also report the min crossing among clause-balanced S
    bal = (6*inside >= m) & (6*outside >= m)
    return m, int(bad.sum()), (int(cross[bal].min()) if bal.any() else None)
for n, D in [(18, 2), (20, 2), (22, 2), (22, 3), (20, 4)]:
    res = [run(n, D, s) for s in range(6)]
    print("n", n, "Delta", D, "m", res[0][0], "m/20 =", res[0][0]/20, "#bad S per seed:", [r[1] for r in res], "min crossing among clause-balanced S:", [r[2] for r in res])
