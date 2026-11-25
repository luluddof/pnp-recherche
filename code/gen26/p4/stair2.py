# combined bound: rank >= (N^2 - |D'|) + max_beta #{(x,y) in D' : pos(y) < t_beta(x)},  D' = dom ∩ im of Theta,
# for the 4 variants (E,C), (E^T,C^T), (C,E), (C^T,E^T). Random validation against exact rank.
import sys, itertools, numpy as np, mc, flint
def Dprime(E, C):
    N = E[0].shape[0]; D = np.zeros((N, N), bool)
    for a in range(N):
        for b in range(N):
            dom = any(np.abs(E[l][:, b]).sum() > 0 and np.abs(C[l][:, a]).sum() > 0 for l in range(len(E)))
            im = any(np.abs(E[l][a, :]).sum() > 0 and np.abs(C[l][b, :]).sum() > 0 for l in range(len(E)))
            D[a, b] = dom and im
    return D
def bound1(E, C):
    N = E[0].shape[0]; Z = (sum(np.abs(m) for m in E) > 0); D = Dprime(E, C); best = 0
    for beta in itertools.permutations(range(N)):
        s = 0
        for x in range(N):
            t = 0
            while t < N and not Z[x, beta[t]]: t += 1
            s += sum(D[x, beta[k]] for k in range(t))
        best = max(best, s)
    return N * N - D.sum() + best
def bounds(E, C):
    T = lambda L: [m.T for m in L]
    return [bound1(E, C), bound1(T(E), T(C)), bound1(C, E), bound1(T(C), T(E))]
if __name__ == "__main__":
    rng = np.random.default_rng(int(sys.argv[1])); viol = 0; gain = 0
    import stair
    for trial in range(int(sys.argv[2])):
        N = int(rng.integers(3, 6)); q = int(rng.integers(1, 5))
        def rpp():
            k = int(rng.integers(0, N + 1)); rows = rng.choice(N, k, replace=False); cols = rng.choice(N, k, replace=False)
            M = np.zeros((N, N), int); M[rows, cols] = rng.choice([-1, 1], k); return M
        E = [rpp() for _ in range(q)]; C = [rpp() for _ in range(q)]
        rk = flint.fmpz_mat((mc.swap(N) - mc.Kmat(N, E, C)).tolist()).rank()
        b = bounds(E, C); b0 = stair.bounds(E, C)
        if max(b) > rk: viol += 1; print("VIOLATION", N, q, rk, b)
        if max(b) > max(b0): gain += 1
    print("trials", sys.argv[2], "violations", viol, "strictly better than plain staircase:", gain)
