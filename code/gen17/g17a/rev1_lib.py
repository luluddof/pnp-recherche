#!/usr/bin/env python3
"""G17a independent checks (reviewer code, not derived from the prover's checks.py)."""
from fractions import Fraction as Fr
import itertools, math, random
INF = 10**9
def v2(x):
    x = Fr(x)
    if x == 0: return INF
    n, d, k = abs(x.numerator), x.denominator, 0
    while n % 2 == 0: n //= 2; k += 1
    while d % 2 == 0: d //= 2; k -= 1
    return k
# ---- polynomials as coefficient lists (low degree first)
def pmul(A, B):
    R = [Fr(0)] * (len(A) + len(B) - 1)
    for i, x in enumerate(A):
        if x == 0: continue
        for j, y in enumerate(B): R[i + j] += x * y
    return R
def padd(A, B):
    n = max(len(A), len(B)); return [(A[i] if i < len(A) else 0) + (B[i] if i < len(B) else 0) for i in range(n)]
def pscale(A, c): return [c * x for x in A]
def ppow(A, k):
    R = [Fr(1)]
    for _ in range(k): R = pmul(R, A)
    return R
def peval(A, x):
    r = Fr(0)
    for c in reversed(A): r = r * x + c
    return r
def pcompose_affine(A, t, al):
    """A(t + al*z) as polynomial in z"""
    R = [Fr(0)]; lin = [Fr(t), Fr(al)]
    for c in reversed(A): R = padd(pmul(R, lin), [c])
    return R
def gauss(A): return min(v2(c) for c in A)
def deg(A):
    d = len(A) - 1
    while d > 0 and A[d] == 0: d -= 1
    return d

def grid(a):
    L = 2 ** a; return L, list(range(-L, L))

