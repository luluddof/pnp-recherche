# Native site-symmetric construction for the P_3 Pauli model: M = (x) M1 (per site 3-cycle on F2^2 \ 0),
# M = w^{D}, D = n_w - n_{w^2}.  C' = P(D), P of degree <= w interpolating w^D on a window of w+1 consecutive D.
# Exact error (eigenvalue count); Pauli support + register projections computed from the explicit operator (k<=4).
import numpy as np, sys, itertools, math
from math import comb
from core import Mperm, pc
k = int(sys.argv[1]); n = 4**k; N = 2**k; msk = N - 1
Mm = Mperm(k).astype(complex)
om = np.exp(2j*np.pi/3)
# per-site d = (M1 - M1^2)/(w - w^2)  (eigen +1 on w-eigvec, -1 on w^2, 0 on 1-eigenspace)
M1 = Mperm(1).astype(complex); d1 = (M1 - M1@M1)/(om - om**2)
I1 = np.eye(4)
def site_op(op, i):  # site i: qubits x_i (bit i) and y_i (bit k+i); M-picture index = x | y<<k
    # build via permutation of tensor factors: basis index -> (x_i,y_i) local index = x_i + 2 y_i
    A = np.zeros((n, n), complex)
    for z in range(n):
        xi = (z >> i) & 1; yi = (z >> (k+i)) & 1; loc = xi + 2*yi
        base = z & ~((1 << i) | (1 << (k+i)))
        for loc2 in range(4):
            if op[loc2, loc] != 0:
                z2 = base | ((loc2 & 1) << i) | ((loc2 >> 1) << (k+i))
                A[z2, z] += op[loc2, loc]
    return A
D = sum(site_op(d1, i) for i in range(k))
assert np.allclose(np.linalg.matrix_power(D, 1) @ Mm, Mm @ D)
# distribution of D over basis (eigen-count): per site +1 w.p. 1/4, -1 1/4, 0 1/2 (dimension count)
cnt = {}
for a in range(k+1):
    for b in range(k+1-a):
        c = k - a - b
        cnt[a-b] = cnt.get(a-b, 0) + math.factorial(k)//(math.factorial(a)*math.factorial(b)*math.factorial(c)) * 2**c
def lagr(W):
    def f(x):
        s = 0
        for i in W:
            t = om**(i % 3)
            for j in W:
                if j != i: t *= (x-j)/(i-j)
            s += t
        return s
    return f
def pauli_proj(C):
    # C = sum_{a} X^a diag(f_a) ; Pauli coeff (a,c) = WHT_z of f_a ; register projections
    S = []
    for a in range(n):
        fa = np.array([C[z ^ a, z] for z in range(n)])
        h = fa.copy(); hh = 1
        while hh < n:
            for i0 in range(0, n, 2*hh):
                x = h[i0:i0+hh].copy(); y = h[i0+hh:i0+2*hh].copy(); h[i0:i0+hh] = x+y; h[i0+hh:i0+2*hh] = x-y
            hh *= 2
        for c in np.nonzero(np.abs(h) > 1e-8*n)[0]: S.append((a, int(c)))
    # a = a1 | a2<<k ; c index: phase (-1)^{c.z}, z = x | y<<k -> c = d | c_y<<k (d on x, c on y)
    r1 = set(); r2 = set(); r3 = set()
    for a, c in S:
        a1, a2 = a & msk, a >> k; d, cy = c & msk, c >> k
        r1.add((a1, cy)); r3.add((a2, d)); r2.add((a1 ^ a2, cy ^ d))
    return len(S), len(r1), len(r2), len(r3)
res = []
for w in range(0, 2*k+1):
    for lo in range(-k, k-w+1):
        W = list(range(lo, lo+w+1)); f = lagr(W)
        err = sum(m for Dv, m in cnt.items() if abs(f(Dv) - om**(Dv % 3)) > 1e-9)
        res.append((w, lo, err))
best = {}
for w, lo, err in res:
    if w not in best or err < best[w][1]: best[w] = (lo, err)
for w in sorted(best):
    lo, err = best[w]
    if k <= 4:
        f = lagr(list(range(lo, lo+w+1)))
        evals, V = np.linalg.eig(D)  # D diagonalisable; build P(D) via polynomial evaluation
        # polynomial evaluation by Horner on the matrix (coefficients from fitting)
        xs = np.arange(lo, lo+w+1); coef = np.polyfit(xs, [om**(x % 3) for x in xs], w) if w > 0 else np.array([om**(lo % 3)])
        C = np.zeros((n, n), complex)
        for cf in coef: C = C @ D + cf*np.eye(n)
        r = np.linalg.matrix_rank(Mm - C, tol=1e-7)
        q, a1, a2, a3 = pauli_proj(C)
        ring_up = min(a1*a2, a1*a3, a2*a3); ring_lo = math.sqrt(a1*a2*a3)
        print(f'k={k} w={w} window=[{lo},{lo+w}] err(formula)={err} rank(num)={r} sigma-rank={N*err} Paulis={q} proj={(a1,a2,a3)} ring in [{ring_lo:.0f},{ring_up}]', flush=True)
    else:
        print(f'k={k} w={w} window=[{lo},{lo+w}] err={err} sigma-rank={N*err}')
