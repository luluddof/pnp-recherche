# power check: how often each bound is tight (mutant bound+1 would be violated)
import sys, numpy as np
sys.argv = [sys.argv[0], sys.argv[1], '0']; exec(open('t.py').read().split('viol = 0')[0])
rng = np.random.default_rng(int(sys.argv[1]) + 100); tight = dict(RSUM=0, CORE=0, STAIR=0, COR4=0, COR4app=0); n = 0
for t in range(900):
    N = int(rng.integers(3, 6)); q = int(rng.integers(1, 7)); mode = t % 3
    E = [rpp(N, int(rng.integers(1, N+1))) for _ in range(q)]
    C = [(e.T * rng.choice([-1, 1], (N, N))) for e in E] if mode == 1 else [rpp(N, int(rng.integers(1, N+1))) for _ in range(q)]
    if mode == 2:
        cur = rk(N, E, C)
        for it in range(40):
            l = int(rng.integers(q)); w = int(rng.integers(2)); k = int(rng.integers(0, N+1))
            E2 = [m.copy() for m in E]; C2 = [m.copy() for m in C]; (E2 if w else C2)[l] = rpp(N, k); v2 = rk(N, E2, C2)
            if v2 <= cur: E, C, cur = E2, C2, v2
    r_ = rk(N, E, C); n += 1
    tight['RSUM'] += (N*N - r_ == sum(abs(e).sum()*abs(c).sum() for e, c in zip(E, C)))
    tight['CORE'] += (r_ == N*N - len(Dp(N, E, C)))
    tight['STAIR'] += (r_ == bnds(N, E, C))
    r = max(int(abs(m).sum()) for m in E + C)
    if r >= 2 and N*N - q*r*r + 1 > 0: tight['COR4app'] += 1; tight['COR4'] += (r_ == N*N - q*r*r + 1)
print(n, tight)
