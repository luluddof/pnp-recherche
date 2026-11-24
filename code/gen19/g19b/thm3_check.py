# check of Theorem 3: charpoly((G(x)G')F) = prod (t-mu_i) prod_{i<j} (t^2 - mu_i mu_j), exactly (sympy, integer matrices)
import sympy as sp, random, itertools
random.seed(1)
def swap(N):
    F = sp.zeros(N*N, N*N)
    for x in range(N):
        for y in range(N): F[x*N+y, y*N+x] = 1
    return F
t = sp.symbols('t')
for N in (2,3,4):
    F = swap(N)
    for trial in range(6):
        G = sp.Matrix(N, N, lambda i,j: random.randint(-3,3)); Gp = sp.Matrix(N, N, lambda i,j: random.randint(-3,3))
        if trial == 0: G = sp.eye(N); Gp = sp.eye(N)
        if trial == 1: Gp = G.T  # some structured case
        Th = sp.kronecker_product(G, Gp) * F
        cp = Th.charpoly(t).as_expr()
        # right side as polynomial in the coefficients of charpoly(GG'): compute via resultants-free route: eigen-free check
        M = G*Gp
        # prod_i (t - mu_i) = charpoly(M)(t);  prod_{i<j}(t^2 - mu_i mu_j) = charpoly of the second exterior power of M at t^2
        idx = list(itertools.combinations(range(N), 2))
        L2 = sp.Matrix(len(idx), len(idx), lambda a,b: M[idx[a][0],idx[b][0]]*M[idx[a][1],idx[b][1]] - M[idx[a][0],idx[b][1]]*M[idx[a][1],idx[b][0]])
        s = sp.symbols('s')
        rhs = M.charpoly(t).as_expr() * (L2.charpoly(s).as_expr().subs(s, t**2) if idx else 1)
        assert sp.expand(cp - rhs) == 0, (N, trial)
        rk = (F - sp.kronecker_product(G, Gp)).rank()
        assert rk >= N*(N-1)//2
        print(N, trial, "charpoly identity OK; rank(F - G(x)G') =", rk, ">=", N*(N-1)//2)
