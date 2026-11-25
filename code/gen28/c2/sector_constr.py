# Exact u=0 sector construction (Theorem-G optimum for S_k-term-wise-symmetric factors), qubit sites, N = 2^k:
# K = sum_{nu != nu'} Pi_nu (x) Pi_nu'  +  sum_nu (Pi_nu (x) Pi_nu)(1_{S(x)S} (x) F_{W_nu})
# predicted p = (N^2 - sum_nu f_nu w_nu^2)/2 ; report osr(K).
import numpy as np, itertools, sys
from math import comb
exec(open('band_check.py').read().split('for k in [4,5,6]:')[0])
def swap(N):
    S = np.zeros((N*N, N*N))
    for x in range(N):
        for y in range(N): S[y*N+x, x*N+y] = 1
    return S
def realign(M, N):
    return M.reshape(N,N,N,N).transpose(0,2,1,3).reshape(N*N,N*N)
def rk(M, tol=1e-8):
    s = np.linalg.svd(M, compute_uv=False); return int((s > tol*max(1,s[0])).sum())
for k in [3,4,5]:
    N = 2**k; S = swap(N); projs = spin_sector_projectors(k)
    perms = list(itertools.permutations(range(k)))
    Pops = [perm_op(k, list(pm)) for pm in perms]
    avg = sum(np.kron(P, P.T) for P in Pops)/len(Pops)   # (1/k!) sum pi (x) pi^{-1}
    K = np.zeros((N*N, N*N))
    pred = 0
    for m in projs:
        for mp in projs:
            if m != mp: K += np.kron(projs[m], projs[mp])
        f = comb(k,m) - (comb(k,m-1) if m>0 else 0); w = k-2*m+1
        PP = np.kron(projs[m], projs[m])
        K += f * (S @ avg @ PP)
        pred += f*w*w
    p = rk(S - K); q = rk(realign(K, N))
    eta = np.log(q)/np.log(N); delta = 2 - np.log(p)/np.log(N)
    print(f"k={k} N={N}: p={p} (pred {(N*N-pred)//2}), osr q={q}; N(N-1)/2={N*(N-1)//2}; eta={eta:.3f} delta={delta:.4f}")
