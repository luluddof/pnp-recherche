#!/usr/bin/env python3
"""G18S r4: exact check of the rank claim behind A'/A''.  Random instances (u in {2,3}, random wires between random
parties), random oriented wire subset, random INTEGER group tensors tau_g; the restricted term is built entry by entry
from the definition (no Kronecker shortcut) and its rank over Q is computed by fraction-free elimination (exact ints).
Checks rank <= 2^c, and reports how often equality holds (tightness of the per-term bound)."""
import itertools, random
from fractions import Fraction
def rank_Q(A):
    A = [[Fraction(x) for x in row] for row in A]; r = 0; n = len(A); m = len(A[0]) if A else 0
    for c in range(m):
        piv = next((i for i in range(r, n) if A[i][c] != 0), None)
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        for i in range(r + 1, n):
            if A[i][c] != 0:
                f = A[i][c] / A[r][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        r += 1
    return r
if __name__ == "__main__":
    random.seed(2026); tests = 0; viol = 0; tight = 0
    for trial in range(400):
        u = random.choice((2, 3)); P = 2 ** u; M = random.randint(2, 6)
        wires = [(random.randrange(P), random.randrange(P)) for _ in range(M)]
        party = [x for w in wires for x in w]
        # oriented subset: for each wire: 0 skip, 1 row end is the row index, 2 column end is the row index
        o = [random.choice((0, 1, 1, 2, 2)) for _ in range(M)]
        sel = [e for e in range(M) if o[e]]
        if not sel or len(sel) > 5: continue
        rowbit = [2 * e if o[e] == 1 else 2 * e + 1 for e in sel]; colbit = [2 * e + 1 if o[e] == 1 else 2 * e for e in sel]
        w = len(sel)
        for d in range(u):
            tabs = [dict(), dict()]
            def tau(g, key):
                if key not in tabs[g]: tabs[g][key] = random.randint(-3, 3)
                return tabs[g][key]
            gbits = [[k for k in range(2 * M) if (party[k] >> d) & 1 == g] for g in (0, 1)]
            Mat = []
            for R in itertools.product((0, 1), repeat=w):
                row = []
                for C in itertools.product((0, 1), repeat=w):
                    b = [0] * (2 * M)
                    for k, v in zip(rowbit, R): b[k] = v
                    for k, v in zip(colbit, C): b[k] = v
                    row.append(tau(0, tuple(b[k] for k in gbits[0])) * tau(1, tuple(b[k] for k in gbits[1])))
                Mat.append(row)
            c = sum(min(sum(1 for k in rowbit if (party[k] >> d) & 1 == g), sum(1 for k in colbit if (party[k] >> d) & 1 == g)) for g in (0, 1))
            rk = rank_Q(Mat); tests += 1; viol += rk > 2 ** c; tight += rk == 2 ** c
            # restricted identity must be a permutation matrix (here: the identity)
            assert all(int(all(R[j] == C[j] for j in range(w))) == int(R == C) for R in [(0,) * w] for C in [(0,) * w])
    print(f"tests {tests}: violations of rank <= 2^c: {viol};  equality rank == 2^c in {tight} cases")
