#!/usr/bin/env python3
"""G18e: EXACT (integer) verification of explicit DF decompositions on tiny instances.
A term is (direction, integer table over all end-bit assignments).  The verifier checks, independently of
how the term was built, (a) that it factors as F(bits in group 0) * G(bits in group 1) (all 2x2 minors
through a pivot vanish) and (b) that the terms sum to I entrywise over Z."""
import itertools

def verify(wires, u, terms, name):
    M = len(wires); nb = 2 * M; party = [p for w in wires for p in w]
    allb = list(itertools.product((0, 1), repeat=nb))
    I = {b: int(all(b[2 * e] == b[2 * e + 1] for e in range(M))) for b in allb}
    tot = {b: 0 for b in allb}; ok = True
    for (d, f) in terms:
        g0 = [k for k in range(nb) if (party[k] >> d) & 1 == 0]; g1 = [k for k in range(nb) if (party[k] >> d) & 1 == 1]
        tab = {}
        for b in allb:
            v = f(b); tot[b] += v
            tab[(tuple(b[k] for k in g0), tuple(b[k] for k in g1))] = v
        piv = next(((a, c) for (a, c), v in tab.items() if v != 0), None)
        if piv is None: continue
        a0, c0 = piv; rows = {a for a, _ in tab}; cols = {c for _, c in tab}
        ok &= all(tab[(a, c)] * tab[(a0, c0)] == tab[(a, c0)] * tab[(a0, c)] for a in rows for c in cols)
    eq = all(tot[b] == I[b] for b in allb)
    print(f"{name}: M={M} u={u} terms={len(terms)}  rank-1 across groups for every term: {ok}  sum == I exactly: {eq}")
    return ok and eq

def two_block(wires, u, dA, dB, X, Y, Z):
    """Id = (Id - SWAP_{XY}) [direction dA, n(n-1)/2 terms] + SWAP_{XY} [direction dB, 1 term], times Id_Z (2^|Z| terms)."""
    def left(e, d):   # bit index of the end of wire e lying in group 0 of direction d
        p, q = wires[e]
        assert ((p >> d) & 1) != ((q >> d) & 1)
        return 2 * e if (p >> d) & 1 == 0 else 2 * e + 1
    for e in X: assert left(e, dA) == left(e, dB)
    for e in Y: assert left(e, dA) != left(e, dB)
    assert len(X) == len(Y)
    xl = [left(e, dA) for e in X]; xr = [l ^ 1 for l in xl]; yl = [left(e, dA) for e in Y]; yr = [l ^ 1 for l in yl]
    vals = list(itertools.product((0, 1), repeat=len(X))); terms = []
    for zv in itertools.product((0, 1), repeat=len(Z)):
        def spect(b, zv=zv): return int(all(b[2 * e] == z and b[2 * e + 1] == z for e, z in zip(Z, zv)))
        for i, j in itertools.combinations(vals, 2):
            def f(b, i=i, j=j, spect=spect):
                L = (tuple(b[k] for k in xl), tuple(b[k] for k in yl)); R = (tuple(b[k] for k in xr), tuple(b[k] for k in yr))
                s = lambda t: (t == (i, j)) - (t == (j, i))
                return 2 * s(L) * s(R) * spect(b) if False else s(L) * s(R) * spect(b)
            terms.append((dA, f))
        def fs(b, spect=spect):
            return int(all(b[k] == b[m] for k, m in zip(xl, yr)) and all(b[k] == b[m] for k, m in zip(yl, xr))) * spect(b)
        terms.append((dB, fs))
    return terms

res = []
w = [(0, 3), (1, 2)];                     res.append(verify(w, 2, two_block(w, 2, 0, 1, [0], [1], []), "u=2 a=2 slots{00,01}      "))
w = [(0, 3), (3, 0), (1, 2)];             res.append(verify(w, 2, two_block(w, 2, 0, 1, [0], [2], [1]), "u=2 a=2 slots{00,11,01}   "))
w = [(0, 3), (3, 0), (1, 2), (2, 1)];     res.append(verify(w, 2, two_block(w, 2, 0, 1, [0, 1], [2, 3], []), "u=2 a=2 full packet (M=4) "))
w = [(0, 7), (1, 6), (2, 5)];             res.append(verify(w, 3, two_block(w, 3, 2, 0, [0], [1], [2]), "u=3 a=3 three pairs       "))
w = [(0, 7), (1, 6), (2, 5), (3, 4)];     res.append(verify(w, 3, two_block(w, 3, 2, 0, [0, 2], [1, 3], []), "u=3 a=3 four pairs (M=4)  "))
print("ALL OK" if all(res) else "FAILURE")
