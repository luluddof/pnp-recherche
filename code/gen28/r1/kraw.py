# Independent exact (p,q) of the singlet family via Pauli-coefficient (Krawtchouk) formula:
# E = (II - XX - YY - ZZ)/4, 1-E = (3II + XX + YY + ZZ)/4.  Coef of a Pauli string of support size u in Pi_m:
# 4^-k sum_j (-1)^j C(u,j) C(k-u,m-j) 3^(k-u-m+j).  q = sum_{u: c_u != 0} C(k,u) 3^u ; p = sum_{m: f(m)!=(-1)^m} C(k,m)3^(k-m)
from fractions import Fraction
from math import comb, ceil, log2
def lag(W):
    def f(m):
        s = Fraction(0)
        for i, a in enumerate(W):
            t = Fraction((-1)**a)
            for j, b in enumerate(W):
                if j != i: t *= Fraction(m - b, a - b)
            s += t
        return s
    return f
def pq(k, W):
    f = lag(W); fv = [f(m) for m in range(k+1)]
    p = sum(comb(k, m)*3**(k-m) for m in range(k+1) if fv[m] != (-1)**m)
    q = 0
    for u in range(k+1):
        c = sum(fv[m]*sum((-1)**j*comb(u, j)*comb(k-u, m-j)*3**(k-u-m+j) for j in range(0, min(u, m)+1) if m-j <= k-u)
                for m in range(k+1))
        if c != 0: q += comb(k, u)*3**u
    return p, q
def rho(n):  # Hurwitz-Radon number
    a = 0
    while n % 2 == 0: n //= 2; a += 1
    c, b = divmod(a, 4); return 8*c + 2**b
if __name__ == "__main__":
    import sys
    # 1) self-dual points d=(k-1)/2, W=[0..d]
    print("k  N  selfdual(p,q)  closed-form  HR-tiling N*ceil(N/rho)  N^2/8  best balanced over all windows")
    for k in range(3, int(sys.argv[1]) + 1):
        N = 2**k
        sd = pq(k, list(range(0, (k-1)//2 + 1))) if k % 2 else None
        cf = sum(comb(k, u)*3**u for u in range((k-1)//2 + 1))
        best = None
        if k <= 16:
            for d in range(k+1):
                for a in range(0, k-d+1):
                    p, q = pq(k, list(range(a, a+d+1)))
                    if best is None or max(p, q) < best[0]: best = (max(p, q), p, q, a, d)
        hr = N*ceil(N/rho(N))
        print(k, N, sd, cf, hr, N*N//8, best, flush=True)
