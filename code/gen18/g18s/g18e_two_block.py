import itertools
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
