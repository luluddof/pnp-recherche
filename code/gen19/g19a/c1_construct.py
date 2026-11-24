#!/usr/bin/env python3
"""G19a c1: exact check of the cut-conditioning identity in model DF and of N = 2^{min_i M_i} on Hamming-1 loads.
Instance: u directions, parties = labels in {0,1}^u, wire e = (p,q): row end bit Y_e at party p, column end bit X_e at q.
Direction d: two groups = bit d of the party.  Bits ordered (Y_0,X_0,Y_1,X_1,...)."""
import itertools, random, sys
import numpy as np
from fractions import Fraction

def target(M):
    T = np.zeros((2,) * (2 * M), dtype=np.int64)
    for b in itertools.product((0, 1), repeat=M):
        T[tuple(x for bb in b for x in (bb, bb))] = 1
    return T

def cut_identity(u, wires, d):
    """returns the list of terms (each = (tau_group0, tau_group1) as full tensors) of the direction-d conditioning identity"""
    M = len(wires)
    cut = [e for e, (p, q) in enumerate(wires) if ((p ^ q) >> d) & 1]
    terms = []
    for z in itertools.product((0, 1), repeat=len(cut)):
        zz = dict(zip(cut, z)); taus = []
        for g in (0, 1):
            t = np.ones((2,) * (2 * M), dtype=np.int64)
            for idx in itertools.product((0, 1), repeat=2 * M):
                v = 1
                for e, (p, q) in enumerate(wires):
                    pg, qg = (p >> d) & 1, (q >> d) & 1
                    if e in zz:
                        if pg == g and idx[2 * e] != zz[e]: v = 0
                        if qg == g and idx[2 * e + 1] != zz[e]: v = 0
                    elif pg == g:   # uncut wire lying inside group g
                        if idx[2 * e] != idx[2 * e + 1]: v = 0
                t[idx] = v
            taus.append(t)
        terms.append(taus)
    return terms, cut

def depends_only_on_group(t, u, wires, d, g):
    """check tau_g is constant in every end bit NOT located in group g of direction d"""
    for e, (p, q) in enumerate(wires):
        for k, party in ((2 * e, p), (2 * e + 1, q)):
            if (party >> d) & 1 != g:
                if not np.array_equal(np.take(t, 0, axis=k), np.take(t, 1, axis=k)): return False
    return True

def app_bound(u, wires):
    """A'': max over oriented wire subsets of min_d  w - sum_g min(e1,e2)   (brute force 3^M)"""
    M = len(wires); best = 0; arg = None
    for o in itertools.product((0, 1, 2), repeat=M):
        sel = [(p, q) if oo == 1 else (q, p) for (p, q), oo in zip(wires, o) if oo]
        if not sel: continue
        vals = []
        for d in range(u):
            c = 0
            for g in (0, 1):
                e1 = sum(1 for (a, b) in sel if (a >> d) & 1 == g); e2 = sum(1 for (a, b) in sel if (b >> d) & 1 == g)
                c += min(e1, e2)
            vals.append(len(sel) - c)
        if min(vals) > best: best = min(vals); arg = o
    return best, arg

def rank_Q(A):
    A = [[Fraction(int(x)) for x in row] for row in A]; r = 0; rows = len(A); cols = len(A[0]) if A else 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i][c] != 0), None)
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        for i in range(rows):
            if i != r and A[i][c] != 0:
                f = A[i][c] / A[r][c]; A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        r += 1
        if r == rows: break
    return r

def flatten(T, wires, orient):
    """restrict to the oriented subset (others' end bits = 0); rows = row-index ends, cols = column-index ends"""
    M = len(wires); sel = [e for e in range(M) if orient[e]]
    rowpos = [2 * e if orient[e] == 1 else 2 * e + 1 for e in sel]; colpos = [2 * e + 1 if orient[e] == 1 else 2 * e for e in sel]
    w = len(sel); A = np.zeros((2 ** w, 2 ** w), dtype=np.int64)
    for i, rb in enumerate(itertools.product((0, 1), repeat=w)):
        for j, cb in enumerate(itertools.product((0, 1), repeat=w)):
            idx = [0] * (2 * M)
            for k in range(w): idx[rowpos[k]] = rb[k]; idx[colpos[k]] = cb[k]
            A[i, j] = T[tuple(idx)]
    return A

cases = [("u=2 Hamming-1, M_i=(2,2)", 2, [(0, 1), (2, 3), (0, 2), (1, 3)]),
         ("u=2 Hamming-1, M_i=(3,2)", 2, [(0, 1), (2, 3), (3, 2), (0, 2), (1, 3)]),
         ("u=3 Hamming-1, M_i=(2,2,2) (6 wires of the 3-cube)", 3, [(0, 1), (6, 7), (0, 2), (5, 7), (0, 4), (3, 7)]),
         ("u=3 Hamming-1, M_i=(1,2,2)", 3, [(0, 1), (0, 2), (5, 7), (0, 4), (3, 7)]),
         ("u=2 mixed distances (not Hamming-1): (0,3),(0,1),(1,3)", 2, [(0, 3), (0, 1), (1, 3)]),
         ("u=2 full antipodal packet", 2, [(0, 3), (3, 0), (1, 2), (2, 1)])]
for name, u, wires in cases:
    M = len(wires); I = target(M)
    Mi = [sum(1 for (p, q) in wires if ((p ^ q) >> d) & 1) for d in range(u)]
    ham1 = all(bin(p ^ q).count('1') == 1 for p, q in wires)
    ok_all = True
    for d in range(u):
        terms, cut = cut_identity(u, wires, d)
        S = sum(t0 * t1 for t0, t1 in terms)
        ok = np.array_equal(S, I) and len(terms) == 2 ** Mi[d] and all(depends_only_on_group(t[g], u, wires, d, g) for t in terms for g in (0, 1))
        ok_all &= ok
    lb, arg = app_bound(u, wires)
    # exact rank certificate for the optimal flattening: rank(I restricted) = 2^w ; rank of each restricted term of the best direction
    A = flatten(I, wires, arg); w = sum(1 for o in arg if o)
    dbest = min(range(u), key=lambda d: Mi[d]); terms, _ = cut_identity(u, wires, dbest)
    rk_terms = max(rank_Q(flatten(t0 * t1, wires, arg)) for t0, t1 in terms)
    print(f"{name}: M={M} M_i={Mi} Hamming-1={ham1} | identity with 2^M_i terms exact & group-factored for every direction: {ok_all}"
          f" | UB log2N<={min(Mi)} | A''={lb} (orientation {arg}, w={w}, rank I|restr={rank_Q(A)}=2^{w}, max rank of a restricted term={rk_terms})"
          f" | tight={lb == min(Mi)}")
# random check of the sandwich  M_i/u-ish <= A'' <= min M_i
random.seed(1); viol = 0; tight = 0; tot = 0; minratio = 9
for trial in range(300):
    u = random.choice((2, 3)); M = random.randint(2, 6)
    wires = []
    while len(wires) < M:
        p, q = random.randrange(2 ** u), random.randrange(2 ** u)
        if p != q: wires.append((p, q))
    Mi = [sum(1 for (p, q) in wires if ((p ^ q) >> d) & 1) for d in range(u)]
    if min(Mi) == 0: continue
    lb, _ = app_bound(u, wires); tot += 1
    if lb > min(Mi): viol += 1
    if lb == min(Mi): tight += 1
    minratio = min(minratio, lb / min(Mi) * u)
print(f"random labellings (u in 2,3; M<=6): {tot} instances, A'' > min M_i: {viol}, A'' = min M_i: {tight}, min of u*A''/min_i M_i = {minratio:.3f} (Lemma Q predicts >= 1 up to rounding -u)")
