# independent test of (1) RSUM, (2) CORE, (3) staircase+CORE, (4) corollary on random MC instances (exact Q rank via flint)
import sys, json, glob, itertools, numpy as np, flint
rng = np.random.default_rng(int(sys.argv[1])); T = int(sys.argv[2])
def rk(N, E, C):
    K = sum(np.kron(e, c) for e, c in zip(E, C)); S = np.zeros((N*N, N*N), int)
    for x in range(N):
        for y in range(N): S[x*N+y, y*N+x] = 1
    return flint.fmpz_mat((S - K).tolist()).rank()
def rpp(N, k):
    rows = rng.choice(N, k, replace=False); cols = rng.choice(N, k, replace=False)
    M = np.zeros((N, N), int); M[rows, cols] = rng.choice([-1, 1], k); return M
def Dp(N, E, C):
    # P_l(a,b) = (row of E_l in col b, row of C_l in col a)
    dom = set(); im = set()
    for e, c in zip(E, C):
        ce = [b for b in range(N) if e[:, b].any()]; cc = [a for a in range(N) if c[:, a].any()]
        re = [x for x in range(N) if e[x].any()]; rc = [y for y in range(N) if c[y].any()]
        dom |= {(a, b) for a in cc for b in ce}; im |= {(x, y) for x in re for y in rc}
    return dom & im
def stair(N, E, C):
    D = Dp(N, E, C); Z = sum(np.abs(e) for e in E) > 0; best = 0
    for beta in itertools.permutations(range(N)):
        s = 0
        for x in range(N):
            t = 0
            while t < N and not Z[x, beta[t]]: t += 1
            s += sum((x, beta[k]) in D for k in range(t))
        best = max(best, s)
    return N*N - len(D) + best
def bnds(N, E, C):
    T_ = lambda L: [m.T for m in L]
    return max(stair(N, E, C), stair(N, T_(E), T_(C)), stair(N, C, E), stair(N, T_(C), T_(E)))
def check(N, E, C, tag=""):
    q = len(E); r_ = rk(N, E, C); v = []
    if sum(abs(e).sum()*abs(c).sum() for e, c in zip(E, C)) < N*N - r_: v.append("RSUM")
    if r_ < N*N - len(Dp(N, E, C)): v.append("CORE")
    if N <= 5 and bnds(N, E, C) > r_: v.append("STAIR")
    r = max(int(abs(m).sum()) for m in E + C)
    if r >= 2 and r_ < N*N - q*r*r + 1: v.append("COR4")
    rE = max(int(abs(m).sum()) for m in E)
    coreE = (rE >= 2 and r_ < N*N - q*rE*rE + 1)   # E-only variant (not the claim)
    if v: print("VIOLATION", tag, N, q, r_, v, [m.tolist() for m in E], [m.tolist() for m in C]); 
    return v, coreE, r_
viol = 0; eonly = 0; minslack = {}
for fn in sorted(glob.glob('/home/user/pnp-recherche/code/gen26/p4/mc_N*.json')):
    d = json.load(open(fn)); E = [np.array(m, int) for m in d['E']]; C = [np.array(m, int) for m in d['C']]; N = E[0].shape[0]
    v, ce, r_ = check(N, E, C, fn); viol += bool(v)
    print(fn.split('/')[-1], "p=", d.get('p'), "rankQ", r_, "core-bound", N*N-len(Dp(N,E,C)), "stair", bnds(N,E,C) if N<=5 else '-', "viol", v)
for t in range(T):
    N = int(rng.integers(3, 6)); q = int(rng.integers(1, 7)); mode = t % 3
    if mode == 0: E = [rpp(N, int(rng.integers(0, N+1))) for _ in range(q)]; C = [rpp(N, int(rng.integers(0, N+1))) for _ in range(q)]
    elif mode == 1:  # transposed-support structure, tends to give low rank
        E = [rpp(N, int(rng.integers(1, N+1))) for _ in range(q)]; C = [(e.T * rng.choice([-1, 1], (N, N))) for e in E]
    else:  # hill-climb: minimise rank by local moves
        E = [rpp(N, int(rng.integers(1, N+1))) for _ in range(q)]; C = [rpp(N, int(rng.integers(1, N+1))) for _ in range(q)]
        cur = rk(N, E, C)
        for it in range(40):
            l = int(rng.integers(q)); w = int(rng.integers(2)); k = int(rng.integers(0, N+1))
            E2 = [m.copy() for m in E]; C2 = [m.copy() for m in C]
            (E2 if w else C2)[l] = rpp(N, k); v2 = rk(N, E2, C2)
            if v2 <= cur: E, C, cur = E2, C2, v2
    v, ce, r_ = check(N, E, C, f"t{t}m{mode}"); viol += bool(v); eonly += ce
print("instances", T, "violations", viol, "E-only-variant violations", eonly)
