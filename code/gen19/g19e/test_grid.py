import sys, numpy as np, itertools
sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from ff import brank, rank1, nullspace
from grid import *
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
def swapN(N):
    F = np.zeros((N*N, N*N), dtype=np.int64)
    for a in range(N):
        for b in range(N): F[b*N+a, a*N+b] = 1
    return F
# --- T0: P built two ways agree (permutation vs product of rectangle swaps)
for N, h in [(2,2),(3,2),(2,3),(3,3)]:
    P = build_P(N, h)
    ops = {((j,k),(k,j)): swapN(N) for j in range(h) for k in range(j+1,h)}
    P2 = kron_slots(ops, N, h)
    assert (P == P2).all()
    for p in (2,3,101):
        assert osr(P, N, h, 0, p) == N**(2*(h-1)) and rank1(P, p) == N**(h*(h-1))
print("T0 ok: P = product of rectangle swaps; rank(P) = N^{h(h-1)}, osr_j(P) = M^2")
# --- T1: product point of the hexagon, N = 2:  F = I - Y, rank Y = 1;  A = -(Y x Y x Y), B = P - A
N, h = 2, 3
P = build_P(N, h); I4 = np.eye(4, dtype=np.int64); Y = I4 - swapN(2)
A = kron_slots({((j,k),(k,j)): Y for j in range(h) for k in range(j+1,h)}, N, h)
for p in (2, 3, 101, 10007):
    B = (P + A) % p            # P = prod(I - Y) ; the YYY term has sign (-1)^3
    r = rank1(A, p); bs = [osr(B, N, h, j, p) for j in range(h)]
    print(f"T1 hexagon N=2 product point mod {p}: rank(A)={r}, osr_j(P-A)={bs}, need >= M^2(1-r/Q)={16*(1-r/4)}; violations {[check(N,h,r,b) for b in bs]}")
# identity: rank(P - 1)
for (N,h) in [(2,3),(3,2),(2,2),(4,2)]:
    P = build_P(N,h); Id = np.eye(P.shape[0], dtype=np.int64)
    for p in (2,3,101):
        r = rank1(P-Id, p); print(f"T1b (N,h)=({N},{h}) mod {p}: rank(P-1)={r}, b=1, bound(a)={bounds_a(N,h,1)}, viol {check(N,h,r,1)}")
# --- T2: random / structured B of small osr, r = rank(P - B)
def rand_B(N, h, b, p, mode):
    S = slots(h); L = len(S); M = N**(h-1); R = N**(L-(h-1))
    B = np.zeros((N**L, N**L), dtype=np.int64)
    P = build_P(N, h)
    for _ in range(b):
        if mode == 'rand':
            Lk = rng.integers(0, p, (M, M)); Mk = rng.integers(0, p, (R, R))
        elif mode == 'perm':
            Lk = np.eye(M, dtype=np.int64)[rng.permutation(M)]; Mk = np.eye(R, dtype=np.int64)[rng.permutation(R)]
        elif mode == 'sparse':
            Lk = (rng.random((M, M)) < 1.5/M).astype(np.int64); Mk = (rng.random((R, R)) < 1.5/R).astype(np.int64)
        elif mode == 'frag':   # a realignment-row of P: L = E_{ab}, M = the matching block of P
            a, c = rng.integers(0, M, 2); Lk = np.zeros((M, M), dtype=np.int64); Lk[a, c] = 1
            T = P.reshape(M, R, M, R); Mk = T[a, :, c, :].copy()
        B += np.kron(Lk, Mk)
    return B % p
tot = 0; worst = {}
for (N, h) in [(2,2),(3,2),(4,2),(2,3)]:
    P = build_P(N, h); M = N**(h-1); Q = N**((h-1)*(h-2))
    for p in (2, 3, 5):
        for mode in ('rand','perm','sparse','frag'):
            for b in range(1, min(M*M, 9)+1):
                for t in range(6 if h == 3 else 20):
                    B = rand_B(N, h, b, p, mode)
                    bb = osr(B, N, h, 0, p); r = rank1(P - B, p)
                    v = check(N, h, r, bb); tot += 1
                    key = (N, h, bb); worst[key] = min(worst.get(key, 10**9), r)
                    assert not v, (N, h, p, mode, bb, r, v)
print(f"T2 ok: {tot} decompositions P = (P-B) + B, no violation.  min r seen per (N,h,b):")
for (N,h) in [(2,2),(3,2),(4,2),(2,3)]:
    print("   ", (N,h), [(b, worst[(N,h,b)], bounds_a(N,h,b)) for (n_,h_,b) in sorted(worst) if (n_,h_)==(N,h)], " (b, min r seen, bound a)")
# --- T3: structured A of small rank, b = osr(P - A)
tot = 0
for (N, h) in [(2,2),(3,2),(4,2),(2,3)]:
    P = build_P(N, h); Dm = P.shape[0]; M = N**(h-1); Q = N**((h-1)*(h-2))
    for p in (2, 3, 5):
        for mode in ('rand', 'Pproj', 'Pcols'):
            for r in range(0, min(Dm, 2*Q*M)+1, max(1, Q//2)):
                for t in range(4 if h==3 else 10):
                    if mode == 'rand':
                        A = (rng.integers(0,p,(Dm,r)) @ rng.integers(0,p,(r,Dm))) % p
                    elif mode == 'Pcols':     # P restricted to r random basis inputs
                        A = np.zeros_like(P); c = rng.choice(Dm, r, replace=False); A[:, c] = P[:, c]
                    else:                     # P times a random rank-r matrix
                        A = (P @ (rng.integers(0,p,(Dm,r)) @ rng.integers(0,p,(r,Dm)))) % p
                    rr = rank1(A, p); b = osr(P - A, N, h, 0, p); v = check(N, h, rr, b); tot += 1
                    assert not v, (N, h, p, mode, rr, b, v)
print(f"T3 ok: {tot} decompositions P = A + (P-A), no violation")
