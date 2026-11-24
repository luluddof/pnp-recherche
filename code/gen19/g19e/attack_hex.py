# Adversarial local search over F_2: minimise rank(P - L (x) M') for the hexagon (2x3, N=2), b = 1.  Theorem F says >= 16.
# Also the crossing N=4 with q=2 (bound 2) and N=5 with q=1 (bound 6, N(N-1)/2 = 10).
import sys, time, random, numpy as np
sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from grid import build_P, slots
def rank2(rows):
    rows = list(rows); r = 0
    while rows:
        x = rows.pop()
        if x:
            r += 1; lb = x & -x
            rows = [y ^ x if y & lb else y for y in rows]
    return r
def setup(N, h):
    P = build_P(N, h); S = slots(h); L = len(S); M = N**(h-1); R = N**(L-(h-1)); Dm = N**L
    # index of out/in: site-0 slots are the first h-1 -> index = i0*R + irest
    Prows = [int(''.join(map(str, P[i][::-1])), 2) for i in range(Dm)]
    return Prows, M, R, Dm
def kron_rows(Ls, Ms, M, R):
    # sum_t L_t (x) M_t as bit rows; L_t: list of M row-masks (M bits), M_t: list of R row-masks (R bits)
    rows = [0]*(M*R)
    for Lm, Mm in zip(Ls, Ms):
        for a in range(M):
            la = Lm[a]
            if not la: continue
            for c in range(R):
                mc = Mm[c]; v = 0
                for u in range(M):
                    if (la >> u) & 1: v |= mc << (u*R)
                rows[a*R + c] ^= v
    return rows
def search(N, h, b, tmax, seed, target):
    random.seed(seed); Prows, M, R, Dm = setup(N, h)
    best_all = Dm; t0 = time.time(); starts = 0
    while time.time() - t0 < tmax:
        starts += 1
        if starts == 1 and b == 1:
            Ls = [[1 << i for i in range(M)]]; Ms = [[1 << i for i in range(R)]]
        else:
            Ls = [[random.getrandbits(M) for _ in range(M)] for _ in range(b)]; Ms = [[random.getrandbits(R) for _ in range(R)] for _ in range(b)]
        cur = rank2([x ^ y for x, y in zip(Prows, kron_rows(Ls, Ms, M, R))]); stall = 0
        while stall < 400 and time.time() - t0 < tmax:
            t = random.randrange(b)
            if random.random() < 0.3:
                a, u = random.randrange(M), random.randrange(M); Ls[t][a] ^= 1 << u; und = ('L', t, a, u)
            else:
                a, u = random.randrange(R), random.randrange(R); Ms[t][a] ^= 1 << u; und = ('M', t, a, u)
            new = rank2([x ^ y for x, y in zip(Prows, kron_rows(Ls, Ms, M, R))])
            if new <= cur:
                stall = 0 if new < cur else stall + 1; cur = new
            else:
                stall += 1
                if und[0] == 'L': Ls[und[1]][und[2]] ^= 1 << und[3]
                else: Ms[und[1]][und[2]] ^= 1 << und[3]
        best_all = min(best_all, cur)
    print(f"(N,h)=({N},{h}) b<={b}: best rank found {best_all} over {starts} starts in {tmax}s; proven bound {target}", flush=True)
    assert best_all >= target
search(2, 3, 1, 70, 1, 16)
search(2, 3, 2, 50, 2, 8)
search(4, 2, 2, 50, 3, 2)
search(5, 2, 1, 30, 4, 6)
search(3, 2, 2, 20, 5, 1)
