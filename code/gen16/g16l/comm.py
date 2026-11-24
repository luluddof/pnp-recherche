# Independent check (G16L). Commutative class: J = sum_z g_z(t) X^z, eigenvalue polys lam_w(t), lam_w(0)=(-1)^|w|,
# deg<=d. rank of 01 block at point k (direction c) = #{edges {w,w+e_c}: lam_w(tau_k) != lam_{w+e_c}(tau_k)}.
# Exact min of sum_k rank = N*2^(l-1) - max feasible subsystem of edge equalities. DFS w/ incremental elimination.
import sys, itertools
from fractions import Fraction
sys.setrecursionlimit(10000)

def bound_new(l, m, d):
    N = sum(m); s = 0
    for x in itertools.product([0,1], repeat=l):
        mu = sum(m[c] for c in range(l) if x[c]); s += max(0, N-d-2*mu)
    return s
def bound_old(l, m, d):
    N = sum(m); mm = (N-d-1)//2; s = 0
    if N-d-1 < 0: return 0
    for x in itertools.product([0,1], repeat=l):
        mu = sum(m[c] for c in range(l) if x[c]); s += max(0, mm+1-mu)
    return s

def solve(l, m, d, p, taus=None, limit=None):
    """return exact min cost (or >=limit if limit given and not beaten)"""
    N = sum(m)
    if p:
        inv = lambda a: pow(a, p-2, p); red = lambda a: a % p
    else:
        inv = lambda a: 1/Fraction(a); red = lambda a: Fraction(a)
    if taus is None: taus = list(range(1, N+1))
    if p: assert len(set(t % p for t in taus)) == N and all(t % p for t in taus)
    dirs = [c for c in range(l) for _ in range(m[c])]
    W = list(itertools.product([0,1], repeat=l)); idx = {w:i for i,w in enumerate(W)}
    nv = len(W)*d
    eqs = []
    for k in range(N):
        c = dirs[k]; t = taus[k]
        for w in W:
            if w[c] == 0:
                w2 = tuple(1-b if i == c else b for i,b in enumerate(w))
                row = [red(0)]*(nv+1)
                for i in range(1, d+1):
                    row[idx[w]*d+i-1] = red(t**i); row[idx[w2]*d+i-1] = red(-t**i)
                sgn = (-1)**sum(w)
                row[nv] = red(-2*sgn)
                eqs.append(row)
    n = len(eqs)
    best = [limit if limit is not None else n+1]
    # order equations: interleave? keep natural
    def reduce(row, basis):
        row = row[:]
        for piv, b in basis:
            f = row[piv]
            if f:
                row = [red(x - f*y) for x, y in zip(row, b)]
        return row
    def dfs(i, basis, cost):
        if cost >= best[0]: return
        if i == n:
            best[0] = cost; return
        r = reduce(eqs[i], basis)
        piv = next((j for j in range(nv) if r[j]), None)
        if piv is None:
            if r[nv]: dfs(i+1, basis, cost+1)
            else: dfs(i+1, basis, cost)
            return
        iv = inv(r[piv]); r = [red(x*iv) for x in r]
        dfs(i+1, basis+[(piv, r)], cost)
        dfs(i+1, basis, cost+1)
    dfs(0, [], 0)
    return best[0]

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "small":
        bad = 0; tight = 0; tot = 0
        for p in [0, 7, 11, 13, 101]:
            for l in [1,2,3]:
                for N in range(1, 7 if l < 3 else 6):
                    if p and N >= p: continue
                    for m in itertools.product(range(N+1), repeat=l):
                        if sum(m) != N or list(m) != sorted(m, reverse=True): continue
                        for d in range(0, N+1):
                            if p == 0 and l == 3 and N*d > 8: continue
                            bn = bound_new(l, m, d); bo = bound_old(l, m, d)
                            T = solve(l, m, d, p)
                            tot += 1
                            flag = "" if T >= bn else "  **VIOLATION**"
                            if T < bn: bad += 1
                            if T == bn: tight += 1
                            print(f"p={p} l={l} m={m} d={d}: T={T} new={bn} old={bo}{flag}", flush=True)
        print("TOTAL", tot, "violations", bad, "tight", tight)
    elif mode == "one":
        l = int(sys.argv[2]); d = int(sys.argv[3]); p = int(sys.argv[4])
        lim = int(sys.argv[5]) if len(sys.argv) > 5 else None
        m = tuple([1]*l)
        print(f"one-per-direction l={l} d={d} p={p}: T={solve(l,m,d,p,limit=lim)} (limit {lim}) new={bound_new(l,m,d)} old={bound_old(l,m,d)}", flush=True)
