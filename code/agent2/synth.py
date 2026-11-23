"""Exact circuit synthesis over the full binary basis B2 (Kojevnikov-Kulikov-Yaroslavtsev style SAT encoding).
Circuit: n inputs, s gates, fan-in 2, any of the 16 binary functions; size = #gates.
Normalization: all gates are 'normal' (g(0,0)=0); if f(0..0)=1 we synthesize NOT f (negation is free in B2).
Symmetry breaking: every non-output gate is used; gates nondegenerate (not const, not projection).
"""
import itertools, sys, time
from pysat.solvers import Cadical153
from pysat.formula import IDPool
from threading import Timer

def tt_bits(f, n):  # f: int truth table, bit t = f(t), input j = bit j of t
    return [(f >> t) & 1 for t in range(1 << n)]

def exists_circuit(tt, n, s, timeout=None):
    """tt: list of 0/1 with tt[0]==0. Returns True/False/None(timeout)."""
    N = 1 << n
    pool = IDPool()
    X = lambda i, t: pool.id(('x', i, t))
    S = lambda i, j, k: pool.id(('s', i, j, k))
    F = lambda i, p, q: pool.id(('f', i, p, q))
    cl = []
    gates = range(n, n + s)
    def val(j, t):  # literal for node j value at row t, or constant
        if j < n:
            return ('c', (t >> j) & 1)  # constant
        return X(j, t)
    for i in gates:
        sel = [S(i, j, k) for j in range(i) for k in range(j + 1, i)]
        cl.append(sel)
        for a in range(len(sel)):
            for b in range(a + 1, len(sel)):
                cl.append([-sel[a], -sel[b]])
        # nondegenerate: not 0 (f01=f10=f11=0), not proj p (f01=0,f10=1,f11=1), not proj q (f01=1,f10=0,f11=1)
        cl.append([F(i, 0, 1), F(i, 1, 0), F(i, 1, 1)])
        cl.append([F(i, 0, 1), -F(i, 1, 0), -F(i, 1, 1)])
        cl.append([-F(i, 0, 1), F(i, 1, 0), -F(i, 1, 1)])
        for j in range(i):
            for k in range(j + 1, i):
                for t in range(1, N):
                    for p in (0, 1):
                        for q in (0, 1):
                            base = [-S(i, j, k)]
                            ok = True
                            for node, bit in ((j, p), (k, q)):
                                v = val(node, t)
                                if isinstance(v, tuple):
                                    if v[1] != bit:
                                        ok = False
                                else:
                                    base.append(-v if bit else v)
                            if not ok:
                                continue
                            if p == 0 and q == 0:
                                cl.append(base + [-X(i, t)])
                            else:
                                cl.append(base + [-X(i, t), F(i, p, q)])
                                cl.append(base + [X(i, t), -F(i, p, q)])
    out = n + s - 1
    for t in range(1, N):
        cl.append([X(out, t)] if tt[t] else [-X(out, t)])
    # every non-output gate used
    for i in range(n, n + s - 1):
        cl.append([S(m, j, k) for m in range(i + 1, n + s) for j in range(m) for k in range(j + 1, m) if i in (j, k)])
    # topological-order symmetry breaking: if gate i+1 does not read gate i, then (j_i,k_i) <=lex (j_{i+1},k_{i+1})
    for i in range(n, n + s - 1):
        pi = [(j, k) for j in range(i) for k in range(j + 1, i)]
        pn = [(j, k) for j in range(i + 1) for k in range(j + 1, i + 1) if i not in (j, k)]
        for (j, k) in pi:
            for (a, b) in pn:
                if (a, b) < (j, k):
                    cl.append([-S(i, j, k), -S(i + 1, a, b)])
    solver = Cadical153(bootstrap_with=cl)
    if timeout:
        timer = Timer(timeout, lambda: solver.interrupt())
        timer.start()
        r = solver.solve_limited(expect_interrupt=True)
        timer.cancel()
    else:
        r = solver.solve()
    model = solver.get_model() if r else None
    solver.delete()
    return r, model, pool

def depends(tt, n):
    return [j for j in range(n) if any(tt[t] != tt[t ^ (1 << j)] for t in range(1 << n))]

def circuit_size(f, n, timeout=None, smax=20, verbose=False):
    tt = tt_bits(f, n)
    if tt[0]:
        tt = [1 - b for b in tt]
    d = len(depends(tt, n))
    if d == 0:
        return 0, True
    if d == 1:
        return 0, True  # x or not x
    lo = d - 1  # connected circuit needs >= d-1 gates
    for s in range(lo, smax + 1):
        t0 = time.time()
        r, _, _ = exists_circuit(tt, n, s, timeout)
        if verbose:
            print(f'  s={s}: {r} ({time.time()-t0:.1f}s)', flush=True)
        if r is None:
            return s, False  # size >= s, unknown
        if r:
            return s, True
    return None, False

if __name__ == '__main__':
    n = int(sys.argv[1]); f = int(sys.argv[2], 0)
    print(circuit_size(f, n, verbose=True))
