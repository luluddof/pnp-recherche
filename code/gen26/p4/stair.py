# check the E-side staircase lemma in MC: rank(I - Theta) >= max_beta sum_x t_beta(x), t = leading zero columns (beta order)
# of row x of the union support Z of the E_l ; 4 variants: E rows, E^T rows (cols), C rows, C^T rows.
import sys, json, itertools, numpy as np, mc, flint
def tbound(Z):
    N = Z.shape[0]; best = 0
    for beta in itertools.permutations(range(N)):
        s = 0
        for x in range(N):
            t = 0
            while t < N and Z[x, beta[t]] == 0: t += 1
            s += t
        best = max(best, s)
    return best
def bounds(E, C):
    ZE = (sum(np.abs(m) for m in E) > 0).astype(int); ZC = (sum(np.abs(m) for m in C) > 0).astype(int)
    return [tbound(ZE), tbound(ZE.T), tbound(ZC), tbound(ZC.T)]
if __name__ == "__main__":
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    viol = 0; tot = 0
    for trial in range(400):
        N = int(rng.integers(3, 6)); q = int(rng.integers(1, 5))
        def rpp():
            k = int(rng.integers(0, N + 1)); rows = rng.choice(N, k, replace=False); cols = rng.choice(N, k, replace=False)
            M = np.zeros((N, N), int); M[rows, cols] = rng.choice([-1, 1], k); return M
        E = [rpp() for _ in range(q)]; C = [rpp() for _ in range(q)]
        rk = flint.fmpz_mat((mc.swap(N) - mc.Kmat(N, E, C)).tolist()).rank()
        b = bounds(E, C); tot += 1
        if max(b) > rk: viol += 1; print("VIOLATION", N, q, rk, b)
    print("random trials", tot, "violations", viol)
    for fn, kE, kC in [(a.split(':')[0], a.split(':')[1], a.split(':')[2]) for a in sys.argv[2:]]:
        d = json.load(open(fn)); E = [np.array(m, dtype=float).astype(int) for m in d[kE]]; C = [np.array(m, dtype=float).astype(int) for m in d[kC]]
        N = E[0].shape[0]; rk = flint.fmpz_mat((mc.swap(N) - mc.Kmat(N, E, C)).tolist()).rank()
        print(fn.split('/')[-1], "rank", rk, "staircase bounds", bounds(E, C))
