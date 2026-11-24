# 2 x h grid permutation operator P and the two ranks, exact mod p.  Written for G19e (independent of hexlib.py).
import itertools, numpy as np, sys
sys.path.insert(0, __file__.rsplit('/',1)[0])
from ff import brank, rank1, nullspace
def slots(h):
    return [(j, k) for j in range(h) for k in range(h) if k != j]   # slot (j,k): at site j, partner k
def build_P(N, h):
    """P as array with axes (out slots..., in slots...); in slot (j,k) carries strand j->k, out slot (k,j) receives it."""
    S = slots(h); L = len(S); D = N**L
    P = np.zeros((D, D), dtype=np.int64)
    pos = {s: i for i, s in enumerate(S)}
    for vals in itertools.product(range(N), repeat=L):
        out = [0]*L
        for (j, k), a in zip(S, vals):
            out[pos[(k, j)]] = a
        i_in = 0; i_out = 0
        for a in vals: i_in = i_in*N + a
        for a in out: i_out = i_out*N + a
        P[i_out, i_in] = 1
    return P
def osr(T, N, h, site, p):
    """operator-Schmidt rank of T (D x D) across site | other sites."""
    S = slots(h); L = len(S)
    A = T.reshape([N]*(2*L))
    mine = [i for i, s in enumerate(S) if s[0] == site]
    rest = [i for i in range(L) if i not in mine]
    perm = mine + [L+i for i in mine] + rest + [L+i for i in rest]
    A = A.transpose(perm).reshape(N**(2*len(mine)), -1)
    return rank1(A, p)
def local2(O, N, h, s, t):
    """embed a two-slot operator O[(o_s,o_t),(i_s,i_t)] acting on slots s,t (identity elsewhere is NOT added: returns factor info)"""
    raise NotImplementedError
def kron_slots(ops, N, h):
    """ops: dict {(s,t): O (N^2 x N^2) on slot pair (s,t)} covering all slots; returns the D x D operator."""
    S = slots(h); L = len(S); pos = {s: i for i, s in enumerate(S)}
    T = None; order = []
    for (s, t), O in ops.items():
        T = O if T is None else np.kron(T, O)
        order += [pos[s], pos[t]]
    A = T.reshape([N]*(2*L))
    # axes currently: out in 'order' ordering, then in in 'order' ordering; move to natural ordering
    inv = [order.index(i) for i in range(L)]
    A = A.transpose(inv + [L + i for i in inv])
    return A.reshape(N**L, N**L)
def bounds_a(N, h, b):
    M = N**(h-1); Q = N**((h-1)*(h-2))
    return max(k*Q*(M - k*b) for k in range(M+1))
def check(N, h, r, b):
    """Theorem C/F inequalities for site rank pair (r,b): returns list of violated ones."""
    M = N**(h-1); Q = N**((h-1)*(h-2)); bad = []
    for k in range(M+1):
        if r < k*Q*(M - k*b): bad.append(('a', k))
        if b*Q < k*(M*Q - k*r): bad.append(('b', k))
    return bad
