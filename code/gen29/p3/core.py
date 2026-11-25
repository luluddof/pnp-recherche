# Core for the P_3 Pauli model.  sigma|u1,u2,u3> = |u3,u1,u2> on (C^2)^{(x)k} x3.
# M-picture: basis |x,y>, x,y in F2^k, idx = x | (y<<k).  M|x,y> = |y, x^y>.
# Pauli g(a1,a2,c,d) := X(a1,a2) Z(c,d):  |x,y> -> (-1)^{c.y + d.x} |x^a1, y^a2>.
# sigma-picture image: [X^a1 Z^c] (x) [X^{a1+a2} Z^{c+d}] (x) [X^a2 Z^d]  (register Paulis P, PQ, Q).
import numpy as np
from flint import nmod_mat
P1, P2 = 1000003, 998244353
def popc(a):
    a = np.asarray(a, dtype=np.int64); c = np.zeros_like(a)
    while np.any(a):
        c += a & 1; a = a >> 1
    return c
def pc(a): return bin(a).count('1')
def rank_mod(A, p):
    A = np.asarray(A, dtype=object) % p
    n, m = A.shape
    return nmod_mat(n, m, [int(v) for v in A.ravel()], p).rank()
def rank2(A):
    return rank_mod(A, P1), rank_mod(A, P2)
def Mperm(k):
    n = 4**k; msk = (1 << k) - 1; Mm = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        x, y = i & msk, i >> k
        Mm[y | ((x ^ y) << k), i] = 1
    return Mm
def pauli_M(k, a1, a2, c, d):
    n = 4**k; msk = (1 << k) - 1; P = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        x, y = i & msk, i >> k
        P[(x ^ a1) | ((y ^ a2) << k), i] = (-1)**(pc(c & y) + pc(d & x))
    return P
def fwht(v):
    v = np.array(v, dtype=object).copy(); h = 1; n = len(v)
    while h < n:
        for i in range(0, n, 2*h):
            a = v[i:i+h].copy(); b = v[i+h:i+2*h].copy()
            v[i:i+h] = a + b; v[i+h:i+2*h] = a - b
        h *= 2
    return v
