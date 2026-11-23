"""Exact circuit complexity over the full binary basis B2 (fan-in 2, all 16 ops, NOT free via op choice;
output may be negated / be an input / constant -> size 0). SAT encoding a la Kojevnikov-Kulikov-Yaroslavtsev / Knuth."""
import itertools, sys, json, time
import numpy as np
from pysat.solvers import Cadical153 as Solver

def tt_inputs(n):
    N = 1 << n
    return [sum(((t >> i) & 1) << t for t in range(N)) for i in range(n)]

def size0(f, n):
    N = 1 << n; full = (1 << N) - 1
    cand = {0, full} | set(tt_inputs(n)) | {full ^ x for x in tt_inputs(n)}
    return f in cand

def exists_circuit(f, n, s):
    """Is there a B2 circuit with exactly s gates (last gate = output, possibly negated) computing f?"""
    N = 1 << n
    cnt = [0]
    def new():
        cnt[0] += 1; return cnt[0]
    cl = []
    # x[i][t]: value of gate i on input t, for t>=1 (t=0 all-zero input: normalized gates output 0)
    x = [[new() for t in range(N)] for i in range(s)]
    # op bits g[i][(a,b)] for (a,b) in {(0,1),(1,0),(1,1)}; g(0,0)=0 (normal gates)
    g = [{ab: new() for ab in [(0,1),(1,0),(1,1)]} for i in range(s)]
    sel = [dict() for i in range(s)]
    for i in range(s):
        m = n + i
        pairs = list(itertools.combinations(range(m), 2))
        for (j, k) in pairs:
            sel[i][(j, k)] = new()
        cl.append([sel[i][p] for p in pairs])
        for p, q in itertools.combinations(pairs, 2):
            cl.append([-sel[i][p], -sel[i][q]])
    def val(node, t):  # returns (literal or constant)
        if node < n:
            return ('c', (t >> node) & 1)
        return ('v', x[node - n][t])
    for i in range(s):
        for t in range(N):
            xi = x[i][t]
            for (j, k), sv in sel[i].items():
                vj, vk = val(j, t), val(k, t)
                for a in (0, 1):
                    for b in (0, 1):
                        for c in (0, 1):
                            # if sel and v_j==a and v_k==b then (x_i == c) iff g(a,b)==c
                            base = [-sv]
                            ok = True
                            for (kind, v), want in ((vj, a), (vk, b)):
                                if kind == 'c':
                                    if v != want: ok = False
                                else:
                                    base.append(-v if want else v)
                            if not ok: continue
                            if (a, b) == (0, 0):
                                if c == 1: cl.append(base + [-xi])
                                continue
                            gv = g[i][(a, b)]
                            cl.append(base + [(xi if c else -xi), (-gv if c else gv)])
        # nondegenerate: op not constant/projection (pure redundancy pruning)
        g01, g10, g11 = g[i][(0,1)], g[i][(1,0)], g[i][(1,1)]
        cl.append([g01, g10, g11])          # not const 0
        cl.append([g01, -g10, -g11])        # not projection a
        cl.append([-g01, g10, -g11])        # not projection b
    # every non-output gate must be used
    for i in range(s - 1):
        users = [sv for i2 in range(i + 1, s) for (j, k), sv in sel[i2].items() if j == n + i or k == n + i]
        cl.append(users)
    # output: last gate, possibly negated. normalized output value on t=0 is 0 so negate iff f(0)=1
    neg = f & 1
    for t in range(N):
        bit = ((f >> t) & 1) ^ neg
        if t == 0:
            if bit != 0: return False
            continue
        cl.append([x[s - 1][t] if bit else -x[s - 1][t]])
    # x[i][0] variables unconstrained-but-irrelevant: force 0
    for i in range(s):
        cl.append([-x[i][0]])
    with Solver(bootstrap_with=cl) as S:
        return S.solve()

def complexity(f, n, maxs=12):
    if size0(f, n): return 0
    for s in range(1, maxs + 1):
        if exists_circuit(f, n, s): return s
    return None
