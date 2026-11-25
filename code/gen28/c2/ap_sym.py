# Alternating projections for SWAP = R + C, rank Fl1(R) <= p, rank Fl2(C) <= q (osr), with optional S_k site symmetry.
# N = d^k. Tensor legs: x (k sites), y (k), u (k), v (k);  SWAP[x,y,u,v] = [x=v][y=u]. Fl1 rows (x,y) cols (u,v); Fl2 rows (x,u) cols (y,v).
# usage: python3 ap_sym.py d k sym(0/1) q p1,p2,... iters restarts
import numpy as np, itertools, sys, time
d, k, sym, q = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
plist = [int(a) for a in sys.argv[5].split(',')]; iters = int(sys.argv[6]); restarts = int(sys.argv[7])
N = d**k; rng = np.random.default_rng(int(time.time()) % 1000)
S = np.zeros((N,N,N,N))
for x in range(N):
    for y in range(N): S[x,y,y,x] = 1
def fl1(T): return T.reshape(N*N, N*N)
def un1(M): return M.reshape(N,N,N,N)
def fl2(T): return T.transpose(0,2,1,3).reshape(N*N, N*N)
def un2(M): return M.reshape(N,N,N,N).transpose(0,2,1,3)
perms = list(itertools.permutations(range(k)))
def symm(T):
    if not sym: return T
    big = T.reshape([d]*(4*k)); acc = np.zeros_like(big)
    for pm in perms:
        ax = []
        for leg in range(4): ax += [leg*k + pm[i] for i in range(k)]
        acc += big.transpose(ax)
    return (acc/len(perms)).reshape(N,N,N,N)
def trunc(M, r):
    U, s, Vt = np.linalg.svd(M, full_matrices=False)
    return (U[:, :r]*s[:r]) @ Vt[:r]
for p in plist:
    best = 1e9
    for rs in range(restarts):
        R = symm(un1(trunc(fl1(S + 0.3*rng.normal(size=S.shape)), p)))
        for it in range(iters):
            C = symm(un2(trunc(fl2(S - R), q)))
            R = symm(un1(trunc(fl1(S - C), p)))
            if it % 50 == 0 or it == iters-1:
                res = np.linalg.norm(S - R - C)
                if res < 1e-9: break
        best = min(best, res)
        if best < 1e-9: break
    print(f"d={d} k={k} N={N} sym={sym} q={q} p={p}: best residual {best:.2e} ({'FEASIBLE' if best < 1e-8 else 'not found'})", flush=True)
