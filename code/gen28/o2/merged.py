# Exact (mod 2 primes) row rank / column rank of the merged type-class pieces A_i = sum_{|Q|=i} R^{Q^c} (x) C^{Q} of the
# k-th tensor power of a piece, vs. entropy (subadditive) and multiplicity-free predictions.
import sys, itertools, math, numpy as np, flint
from pieces import *
PR = [1000003, 1000033]
def rk(A):
    out = []
    for pr in PR:
        out.append(flint.nmod_mat(A.shape[0], A.shape[1], [int(v) % pr for v in A.ravel()], pr).rank())
    return out
def kron_list(ms):
    out = ms[0]
    for m in ms[1:]: out = np.kron(out, m)
    return out
def merged(Rm, Cm, k, i):
    A = None
    for Q in itertools.combinations(range(k), i):
        T = kron_list([Cm if s in Q else Rm for s in range(k)])
        A = T if A is None else A + T
    return A
def run(name, N, R, C, ks):
    N2 = N*N
    Rr, Cr, Rc, Cc = rowflat(R), rowflat(C), colflat(R), colflat(C)
    p, rC = rk(Rr)[0], rk(Cr)[0]; cR, q = rk(Rc)[0], rk(Cc)[0]
    print(f"== {name}: N={N} p={p} q={q}  rowrank(C)={rC} colrank(R)={cR}", flush=True)
    for k in ks:
        for i in range(k+1):
            r = rk(merged(Rr, Cr, k, i)); c = rk(merged(Rc, Cc, k, i)); b = math.comb(k, i)
            print(f"  k={k} i={i}: merged row rank {r}  | entropy {b*p**(k-i)*N2**i} (true-subranks {b*p**(k-i)*rC**i})"
                  f" | mult-free {p**(k-i)*N2**i} (true-subranks {p**(k-i)*rC**i})", flush=True)
            print(f"            merged col rank {c}  | entropy {b*q**i*N2**(k-i)} (true-subranks {b*q**i*cR**(k-i)})"
                  f" | mult-free {q**i*N2**(k-i)} (true-subranks {q**i*cR**(k-i)})", flush=True)
if __name__ == '__main__':
    which = sys.argv[1]; ks = [int(a) for a in sys.argv[2].split(',')]
    if which == 'oct': N, (R, C, s) = octonion()
    elif which.startswith('F1_'): N, (R, C, s) = F1(int(which[3:]))
    else: N, (R, C, s) = load_json(which)
    run(which, N, R, C, ks)
