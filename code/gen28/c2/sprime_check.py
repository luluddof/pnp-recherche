# Check Theorem S' (sector/grading bound) on explicit qubit-site constructions (d=2).
# Bound: rank(SWAP-K) >= N^2 * Pr[m2 - m1 > 2u], m = second row of the Schur-Weyl shape (k-m, m) of a uniform word,
# valid when every factor A_l, B_l of K = sum A_l (x) B_l has "S_k-degree" <= u (e.g. u-local Pauli strings).
import numpy as np, itertools, sys
from math import comb
rng = np.random.default_rng(1)
I2 = np.eye(2); X = np.array([[0,1],[1,0]]); Y = np.array([[0,-1j],[1j,0]]); Z = np.diag([1.,-1])
P = [I2, X, Y, Z]

def kron_list(ms):
    out = np.array([[1.0+0j]])
    for m in ms: out = np.kron(out, m)
    return out

def swap(N):
    S = np.zeros((N*N, N*N))
    for x in range(N):
        for y in range(N):
            S[y*N+x, x*N+y] = 1
    return S

def sw_prob(k):
    # Pr[m] for m = 0..k//2 : f^{(k-m,m)} * dim W / 2^k
    pr = {}
    for m in range(k//2+1):
        f = comb(k,m) - (comb(k,m-1) if m>0 else 0)
        pr[m] = f*(k-2*m+1)/2**k
    return pr

def bound(k, u):
    pr = sw_prob(k); N = 2**k
    return N*N*sum(pr[a]*pr[b] for a in pr for b in pr if b - a > 2*u)

def pauli_strings(k, w):
    for supp_size in range(w+1):
        for S in itertools.combinations(range(k), supp_size):
            for labs in itertools.product([1,2,3], repeat=supp_size):
                lab = [0]*k
                for s,l in zip(S,labs): lab[s] = l
                yield tuple(lab)

def rand_local_op(k, w, nterms):
    strs = list(pauli_strings(k, w))
    idx = rng.choice(len(strs), size=min(nterms, len(strs)), replace=False)
    A = 0
    for i in idx:
        A = A + rng.normal()*kron_list([P[a] for a in strs[i]])
    return A

def rank(M, tol=1e-8):
    s = np.linalg.svd(M, compute_uv=False)
    return int((s > tol*max(1,s[0])).sum())

if __name__ == "__main__":
    for k in [2,3,4]:
        N = 2**k; S = swap(N)
        for u in range(0, k):
            b = bound(k, u)
            worst = None
            # (a) random u-local factor products, q terms
            for trial in range(6):
                q = rng.integers(1, 6)
                K = sum(np.kron(rand_local_op(k,u,8), rand_local_op(k,u,8)) for _ in range(q))
                # rescale K to make it as close to SWAP as possible (least squares on scalar) -- harmless
                p = rank(S - K)
                worst = p if worst is None else min(worst, p)
            # (b) Pauli truncation of SWAP: K = (1/N) sum_{|P|<=u} P (x) P^dag   (and T3 window polynomials below)
            K = 0
            for lab in pauli_strings(k, u):
                Pm = kron_list([P[a] for a in lab]); K = K + np.kron(Pm, Pm.conj())/N
            ptr = rank(S - K)
            # (c) best symmetric window family f(M): eigenvalue on m-singlet Bell states; Pauli weight <= deg f = u
            bestw = None
            for lo in range(0, k-u+1):
                W = range(lo, lo+u+1)
                p = sum(comb(k,m)*3**(k-m) for m in range(k+1) if m not in W)
                bestw = p if bestw is None else min(bestw, p)
            print(f"k={k} N={N} u={u}: bound={b:.1f}  random-local min p={worst}  plain-trunc p={ptr}  best window family p={bestw}")
            sys.stdout.flush()
