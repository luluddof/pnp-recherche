# General (non-commutative) class, l=2, over GF(p): EXACT min over all J of sum_k rank(01 block of J(tau_k)).
# dir-0 block [[A1,S],[T1,A2]], dir-1 block [[C1,S],[T2,C2]]; S,T* have constant term 1, A*,C* constant 0, deg<=d.
import sys, itertools, numpy as np
def bound_new(m, d):
    N = sum(m); s = 0
    for x in itertools.product([0,1], repeat=2):
        mu = sum(m[c] for c in range(2) if x[c]); s += max(0, N-d-2*mu)
    return s
def polys(p, d, c0):
    # all polys with const c0 and deg<=d : array (p^d, d+1)
    cs = np.array(list(itertools.product(range(p), repeat=d)), dtype=np.int64).reshape(p**d, d)
    return np.concatenate([np.full((len(cs),1), c0, dtype=np.int64), cs], axis=1)
def evals(P, taus, p):
    # (npoly, npts)
    V = np.array([[pow(t, i, p) for t in taus] for i in range(P.shape[1])], dtype=np.int64)
    return (P @ V) % p
def run(p, d, m):
    N = sum(m); taus = list(range(1, N+1)); assert N < p
    t0 = taus[:m[0]]; t1 = taus[m[0]:]
    P0 = polys(p, d, 0); P1 = polys(p, d, 1)
    def f(pts):
        # returns array over S-polys of min cost
        if not pts: return np.zeros(len(P1), dtype=np.int64)
        A = evals(P0, pts, p); T = evals(P1, pts, p); S = evals(P1, pts, p)
        n = len(P0)
        out = np.empty(len(P1), dtype=np.int64)
        # precompute A1*A2 over pairs: (n,n,npts)
        AA = (A[:, None, :] * A[None, :, :]) % p
        nzA = (A[:, None, :] != 0) | (A[None, :, :] != 0)
        for si in range(len(P1)):
            best = 10**9
            for ti in range(len(P1)):
                ST = (S[si] * T[ti]) % p
                det = (AA - ST[None, None, :]) % p
                anynz = nzA | ((S[si] != 0) | (T[ti] != 0))[None, None, :]
                rk = np.where(det != 0, 2, np.where(anynz, 1, 0)).sum(axis=2)
                b = rk.min()
                if b < best: best = b
            out[si] = best
        return out
    tot = f(t0) + f(t1)
    return int(tot.min())
if __name__ == "__main__":
    bad = 0; tot = 0; tight = 0
    cases = []
    for p, dmax, Nmax in [(3,3,2),(5,2,4),(7,2,6),(5,3,4),(11,1,8),(13,2,5)]:
        for N in range(1, Nmax+1):
            for m0 in range(N, -1, -1):
                m = (m0, N-m0)
                if m[0] < m[1]: continue
                for d in range(0, dmax+1):
                    if p**(3*d) > 3_000_000: continue
                    cases.append((p,d,m))
    seen = set()
    for (p,d,m) in cases:
        if (p,d,m) in seen: continue
        seen.add((p,d,m))
        T = run(p, d, m); bn = bound_new(m, d); tot += 1
        if T < bn: bad += 1
        if T == bn: tight += 1
        print(f"GF({p}) l=2 m={m} d={d}: exact min sum rank B01 = {T}, bound = {bn}" + ("  **VIOLATION**" if T < bn else ""), flush=True)
    print("TOTAL", tot, "violations", bad, "tight", tight)
