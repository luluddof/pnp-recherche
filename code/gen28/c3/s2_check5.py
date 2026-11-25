# Theorem S'' check: K = g(T) (function of antisym-site pattern), n=2, k<=5.
# Explicit: tau_i = SWAP_i - 1/2 on 2k qubits; K = sum_U c_U prod tau_i ; compare rank(SWAP-K), osr(K) with pattern formulas.
import numpy as np, itertools, sys
from modrank import rank_mod, solve_mod
PR = [2147483629, 1000000007]
def site_swap_perm(k, i):
    N = 2**k; idx = np.arange(N*N); x = idx // N; y = idx % N
    bx = (x >> i) & 1; by = (y >> i) & 1
    x2 = x ^ ((bx ^ by) << i); y2 = y ^ ((bx ^ by) << i)
    return x2 * N + y2
def run(k, trials, rng):
    N = 2**k; D = N*N
    perms = [site_swap_perm(k, i) for i in range(k)]
    full = np.arange(D); xx = full // N; yy = full % N; fullswap = yy * N + xx
    subsets = [U for r in range(k+1) for U in itertools.combinations(range(k), r)]
    pats = list(itertools.product([0,1], repeat=k))
    w = {T: 3**(k-sum(T)) for T in pats}
    for t in range(trials):
        # random support family, random agreement set Z, solve for c mod p (same rational solution for both primes? use integer RHS)
        supp = [U for U in subsets if len(U)<=2] if t==0 else [U for U in subsets if len(U)<=2 and rng.random()<0.8]+[U for U in subsets if len(U)==3 and rng.random()<0.3]
        if not supp: continue
        Z = [T for T in pats if sum(T)<=2][:len(supp)] if t==0 else [T for T in sorted(pats,key=lambda T:(sum(T),rng.random()))][:len(supp)-1]
        res = []
        for p in PR:
            inv2 = pow(2, p-2, p)
            tv = {0: inv2 % p, 1: (-3*inv2) % p}   # tau on sym / anti
            def mono(U, T):
                v = 1
                for i in U: v = v * tv[T[i]] % p
                return v
            if Z:
                M = np.array([[mono(U, T) for U in supp] for T in Z], dtype=np.int64)
                b = np.array([(-1)**sum(T) % p for T in Z], dtype=np.int64)
                x0 = solve_mod(M, b, p)
                if x0 is None: res.append(None); continue
                # add random nullspace component: pick random c satisfying system: solve with random extra pins
                c = x0
            else:
                c = np.array([rng.integers(1, 1000) for _ in supp])
            # pattern formulas
            g = {T: sum(int(c[j]) * mono(U, T) for j, U in enumerate(supp)) % p for T in pats}
            pform = sum(w[T] for T in pats if g[T] != (-1)**sum(T) % p)
            qform = sum(3**len(U) for j, U in enumerate(supp) if c[j] % p)
            # explicit
            K = np.zeros((D, D), dtype=np.int64)
            for j, U in enumerate(supp):
                if c[j] % p == 0: continue
                Mx = np.eye(D, dtype=np.int64)
                for i in U:
                    S = np.zeros((D, D), dtype=np.int64); S[perms[i], np.arange(D)] = 1
                    T_ = (S - inv2 * np.eye(D, dtype=np.int64)) % p
                    Mx = (Mx @ T_) % p
                K = (K + int(c[j]) * Mx) % p
            SW = np.zeros((D, D), dtype=np.int64); SW[fullswap, np.arange(D)] = 1
            prank = rank_mod((SW - K) % p, p)
            # realignment: K[(u,v),(x,y)] -> R[(u,x),(v,y)]
            K4 = K.reshape(N, N, N, N)  # u v x y
            R = K4.transpose(0, 2, 1, 3).reshape(D, D)
            qrank = rank_mod(R, p)
            res.append((pform, qform, prank, qrank))
        print(k, t, len(supp), len(Z), res, flush=True)

rng=np.random.default_rng(5)
run(5,3,rng)
