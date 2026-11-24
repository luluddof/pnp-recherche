#!/usr/bin/env python3
"""G18a run1: sanity against G17a numbers, then single-packet census."""
import sys, time, random
from g18lib import *

def analyse(a, packets, us, sample=None, do_bb=True, balanced=True, out=sys.stdout, cache={}, stable=400, do_sat=True):
    t0 = time.time(); st = cache.setdefault(('st', a), Setup(a)); L = st.L
    M = sum(len(p) for p in packets); N = M * (2 * L - 1)
    basis = None; n = 1; exh = True
    for p in packets:
        key = (a, tuple(p), sample)
        if key not in cache:
            cache[key] = packet_lattice(st, p, sample=sample, rng=random.Random(1), log=sys.stderr, stable=stable)
        b, ng, ex = cache[key]; exh &= ex; npk = len(p) * (2 * L - 1) + 1
        if basis is None: basis, n = b, npk
        else: basis = product_lattice(basis, n, b, npk); n = n + npk - 1
    mats = cache.setdefault(('mats', a, N), taylor_mats(st, N))
    Ph = phi(st, basis, N, M * st.e, mats)
    full = np.concatenate([Ph[t] for t in st.S], axis=1)
    if do_sat: sat, vals = saturate(full); r = len(vals); idx = sum(vals)
    else: sat = full; r = len(basis); idx = -1
    res = dict(a=a, packets=packets, M=M, N=N, rank=r, sat_index_log2=idx, exhaustive=exh)
    if do_bb:
        if ('bb', a, N) not in cache: cache[('bb', a, N)] = bb_lattice(st, N, mats)
        BB, bbloss = cache[('bb', a, N)]
    for u in us:
        I_st = Instance(st, full, N, u); I_sat = Instance(st, sat, N, u)
        row = dict(u=u, formula=formula_D9(a, u, M), wit=witness_g(a, u, M))
        row['ST_unif'] = I_st.need_additive(); row['SAT_unif'] = I_sat.need_additive()
        if do_bb:
            I_bb = Instance(st, BB, N, u); row['BB_unif'] = I_bb.need_additive()
        if balanced and u > 1:
            row['ST_bal'] = I_st.need_balanced(row['ST_unif'])
            if do_bb: row['BB_bal'] = I_bb.need_balanced(row['BB_unif'])
        res[u] = row
        print(f"a={a} L={L} packets={packets} M={M} N={N} rank={r}/{N+1} log2[sat:ST]={idx} exh={exh} u={u}: "
              f"formula(loss_i, exact D.9 need, paper bound)={row['formula']} witness g={row['wit']} | "
              f"uniform need: BB={row.get('BB_unif')} sat(ST)={row['SAT_unif']} ST={row['ST_unif']}"
              + (f" | balanced m (p_i=m*2^(u-1-i)): BB={row.get('BB_bal')} ST={row['ST_bal']}" if 'ST_bal' in row else "")
              + f"  [{time.time()-t0:.1f}s]", file=out, flush=True)
    return res

if __name__ == "__main__":
    which = sys.argv[1]
    if which == "sanity":
        analyse(1, [[0, 1]], [1])
        analyse(1, [[0, 1], [0, 1]], [1])
        analyse(2, [[0, 1]], [1, 2])
        analyse(2, [[0, 1, 2]], [1, 2])
        analyse(2, [[0, 1, 3]], [1, 2])
        analyse(2, [[0, 1, 2, 3]], [1, 2])
    if which == "a2single":
        for sl in ([0], [0, 1], [0, 2], [0, 3], [1, 2], [0, 1, 2], [0, 1, 3], [0, 2, 3], [0, 1, 2, 3]):
            analyse(2, [sl], [1, 2])
    if which == "a1multi":
        for k in range(1, 6): analyse(1, [[0, 1]] * k, [1])
        for k in range(1, 6): analyse(1, [[0]] * k, [1])
        analyse(1, [[0, 1], [0]], [1]); analyse(1, [[0, 1], [0, 1], [0]], [1]); analyse(1, [[0, 1], [0], [0]], [1])
