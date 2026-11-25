# Theorem S'' formula check, qudit sites (n,k), explicit K = sum_U c_U prod_{i in U}(SWAP_i - 1/n) via permutation expansion.
import numpy as np, itertools, sys
from modrank import rank_mod, solve_mod
PR = [2147483629, 1000000007]
def run(n, k, trials, seed):
    rng = np.random.default_rng(seed)
    N = n**k; D = N*N
    idx = np.arange(D); x = idx // N; y = idx % N
    dig = lambda z, i: (z // n**i) % n
    def swapW(W):  # permutation: swap digits i in W between x and y
        x2 = x.copy(); y2 = y.copy()
        for i in W:
            dx = dig(x, i); dy = dig(y, i)
            x2 = x2 + (dy - dx) * n**i; y2 = y2 + (dx - dy) * n**i
        return x2 * N + y2
    subsets = [U for r in range(k+1) for U in itertools.combinations(range(k), r)]
    pats = list(itertools.product([0,1], repeat=k))
    ws, wa = n*(n+1)//2, n*(n-1)//2
    w = {T: wa**sum(T) * ws**(k-sum(T)) for T in pats}
    perms = {W: swapW(W) for W in subsets}
    for t in range(trials):
        dmax = rng.integers(1, k+1)
        supp = [U for U in subsets if len(U) <= dmax and rng.random() < 0.85]
        if len(supp) < 2: continue
        Z = sorted(pats, key=lambda T: (sum(T), rng.random()))[:len(supp)-1]
        res = []
        for p in PR:
            invn = pow(n, p-2, p)
            tv = {0: (1 - invn) % p, 1: (-1 - invn) % p}
            mono = lambda U, T: int(np.prod([tv[T[i]] for i in U], dtype=object)) % p if U else 1
            M = np.array([[mono(U, T) for U in supp] for T in Z], dtype=np.int64)
            b = np.array([(-1)**sum(T) % p for T in Z], dtype=np.int64)
            c = solve_mod(M, b, p)
            if c is None: res.append(None); continue
            g = {T: sum(int(c[j]) * mono(U, T) for j, U in enumerate(supp)) % p for T in pats}
            pform = sum(w[T] for T in pats if g[T] != (-1)**sum(T) % p)
            qform = sum((n*n-1)**len(U) for j, U in enumerate(supp) if c[j] % p)
            # explicit: prod_{i in U}(S_i - invn) = sum_{W subset U} (-invn)^{|U-W|} S_W
            aW = {}
            for j, U in enumerate(supp):
                if c[j] % p == 0: continue
                for r in range(len(U)+1):
                    for W in itertools.combinations(U, r):
                        aW[W] = (aW.get(W, 0) + int(c[j]) * pow(-invn % p, len(U) - r, p)) % p
            K = np.zeros((D, D), dtype=np.int64)
            for W, a in aW.items():
                K[perms[W], idx] = (K[perms[W], idx] + a) % p
            SW = np.zeros((D, D), dtype=np.int64); SW[perms[tuple(range(k))], idx] = 1
            prank = rank_mod((SW - K) % p, p)
            R = K.reshape(N, N, N, N).transpose(0, 2, 1, 3).reshape(D, D)
            qrank = rank_mod(R, p)
            res.append((pform, qform, prank, qrank))
        ok = all(r is not None and r[0] == r[2] and r[1] == r[3] for r in res)
        print(f"n={n} k={k} t={t} |supp|={len(supp)} {res} {'OK' if ok else 'MISMATCH'}", flush=True)
for (n, k, tr) in [(2, 5, 4), (3, 2, 5), (3, 3, 3)]:
    run(n, k, tr, seed=n*10+k)
