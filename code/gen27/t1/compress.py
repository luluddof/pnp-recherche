# Row compression: keep a set of rows on which [all columns | rhs] has full rank r_all (projection injective on the span).
import sys, time, pickle; sys.dont_write_bytecode = True
import numpy as np, flint
def compress(data, Dmax=None):
    P = data['P']; t = time.time()
    keys = [k for k in data['cols'] if Dmax is None or sum(k[0]) <= Dmax]
    A = np.stack([data['cols'][k] for k in keys] + data['rhs'], axis=1) % P     # m x (ncols + n)
    m, nc = A.shape
    MT = flint.nmod_mat(nc, m, A.T.ravel().tolist(), P)
    R, r = MT.rref()
    piv = []; i = 0
    # pivot columns of rref(M^T) = independent rows of A
    Rn = np.array([[int(R[i2, j]) for j in range(m)] for i2 in range(r)], dtype=np.int64) if r else np.zeros((0, m))
    for i2 in range(r):
        nz = np.nonzero(Rn[i2])[0]; piv.append(int(nz[0]))
    rows = np.array(piv)
    out = dict(data); out['cols'] = {k: data['cols'][k][rows] for k in keys}; out['rhs'] = [v[rows] for v in data['rhs']]
    out['rows_kept'] = rows; out['r_all'] = r
    print(f"compress: m={m} ncols={nc} r_all={r}  ({time.time()-t:.0f}s)", flush=True)
    return out
if __name__ == "__main__":
    d = pickle.load(open(sys.argv[1], 'rb')); Dmax = int(sys.argv[2]) if len(sys.argv) > 2 else None
    out = compress(d, Dmax); fn = sys.argv[1].replace('.pkl', f'_cmp{"" if Dmax is None else Dmax}.pkl')
    pickle.dump(out, open(fn, 'wb')); print(fn)
