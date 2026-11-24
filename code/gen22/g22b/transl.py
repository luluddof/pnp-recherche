# Translation-invariant ansatz for ST arrays on N = 2^n: cell (a, a^c) gets colour (B(c), a ^ f(c)); blocks B: F_2^n -> [nb],
# blocks 0..nG-1 are S (G) palette, the rest T (H).  Signs: solved afterwards as an F_2 linear system.  Exact verify by shape.verify_shape.
import sys, itertools, time, json
import numpy as np
from pysat.solvers import Cadical153
from pysat.card import CardEnc
from shape import verify_shape, show
from convert import intercal_graph
def solve(n, nG, nH, size=None, maxtries=50):
    N = 2**n; nb = nG + nH; vid = itertools.count(1)
    F = [[next(vid) for v in range(N)] for c in range(N)]
    B = [[next(vid) for i in range(nb)] for c in range(N)]
    cl = []
    for c in range(N):
        cl.append(F[c][:]); cl.append(B[c][:])
        for v, w in itertools.combinations(range(N), 2): cl.append([-F[c][v], -F[c][w]])
        for i, j in itertools.combinations(range(nb), 2): cl.append([-B[c][i], -B[c][j]])
    top = next(vid)
    if size:
        for i in range(nb):
            enc = CardEnc.equals(lits=[B[c][i] for c in range(N)], bound=size, top_id=top); top = max(top, enc.nv) + 1; cl += enc.clauses
    vid = itertools.count(top + 1)
    SB = {}
    for c1, c2 in itertools.combinations(range(N), 2):
        s = next(vid); SB[(c1, c2)] = s   # s <-> same block
        for i in range(nb):
            cl.append([-B[c1][i], -B[c2][i], s]); cl.append([-s, -B[c1][i], B[c2][i]])
    def sb(x, y): return SB[(min(x, y), max(x, y))]
    for c1, c2 in itertools.combinations(range(N), 2):
        s = sb(c1, c2)
        for v in range(N):
            cl.append([-s, -F[c1][v], -F[c2][v]])                       # f injective on blocks
            w = v ^ c1 ^ c2
            cl.append([-s, -F[c1][v], -F[c2][w]])                       # c -> c ^ f(c) injective on blocks
        for v1 in range(N):
            for v2 in range(N):
                d = v1 ^ v2
                if d == 0: continue
                c, cp = d ^ c2, d ^ c1
                pre = [-s, -F[c1][v1], -F[c2][v2]]
                cl.append(pre + [sb(c, cp)])
                for i in range(nb):     # opposite palette: block of c1 (i) S  => block of c T, etc.
                    for j in range(nb):
                        if (i < nG) == (j < nG): cl.append(pre + [-B[c1][i], -B[c][j]])
                for u in range(N):      # f(c) ^ f(cp) = d
                    cl.append(pre + [-F[c][u], F[cp][u ^ d]])
    cl.append([B[0][0]])
    S = Cadical153(bootstrap_with=cl)
    for t in range(maxtries):
        if not S.solve(): return None, t
        m = set(l for l in S.get_model() if l > 0)
        f = [next(v for v in range(N) if F[c][v] in m) for c in range(N)]
        bl = [next(i for i in range(nb) if B[c][i] in m) for c in range(N)]
        cells = sign_solve(N, f, bl, nG)
        if cells: return cells, t
        S.add_clause([-F[c][f[c]] for c in range(N)] + [-B[c][bl[c]] for c in range(N)])
    return None, maxtries
def sign_solve(N, f, bl, nG):
    col = [[None]*N for _ in range(N)]
    for a in range(N):
        for b in range(N):
            c = a ^ b; col[a][b] = ('S' if bl[c] < nG else 'T', bl[c] * N + (a ^ f[c]))
    rows = []
    for a, a2 in itertools.combinations(range(N), 2):
        for b in range(N):
            for b2 in range(N):
                if b != b2 and col[a][b] == col[a2][b2]:
                    assert col[a][b2] == col[a2][b] and col[a][b2][0] != col[a][b][0]
                    v = 0
                    for (x, y) in ((a, b), (a2, b2), (a, b2), (a2, b)): v ^= 1 << (x*N + y)
                    rows.append((v, 1))
    # gaussian elimination over F_2 (python ints as bitsets)
    piv = {}
    for v, r in rows:
        while v:
            h = v.bit_length() - 1
            if h in piv: pv, pr = piv[h]; v ^= pv; r ^= pr
            else: piv[h] = (v, r); break
        if v == 0 and r == 1: return None
    x = 0
    for h in sorted(piv):
        pv, pr = piv[h]; val = pr ^ (bin(pv & x & ~(1 << h)).count('1') & 1)
        if val: x |= 1 << h
    return [[(col[a][b][0], col[a][b][1], -1 if (x >> (a*N + b)) & 1 else 1) for b in range(N)] for a in range(N)]
if __name__ == "__main__":
    n, nG, nH = map(int, sys.argv[1:4]); size = int(sys.argv[4]) if len(sys.argv) > 4 else None
    t0 = time.time(); cells, t = solve(n, nG, nH, size)
    N = 2**n
    print(f"N={N} translation ansatz blocks G={nG} H={nH} size={size}: {'FOUND' if cells else 'none'} after {t} colourings, {time.time()-t0:.1f}s", flush=True)
    if cells:
        if N <= 16: show(N, cells)
        ok, p, q = verify_shape(N, cells)
        json.dump({"N": N, "p": p, "q": q, "cells": cells}, open(f"transl_{N}_{p}_{q}.json", "w"))
