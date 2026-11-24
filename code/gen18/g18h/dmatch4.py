"""D=4, lambda=1: is T_eps (8 antipodal eps-pairs on the 16 corners of Q_4) in the span of the eps-products over
perfect matchings of Q_4 by edges (272 of them; each such term sits exactly at the floor 1/4)?  Exact (flint, over Q)."""
import numpy as np, flint, sys
from dmatch import matchings, tensor
D = 4; n = 16; full = 15
ham = lambda a, b: bin(a ^ b).count("1")
anti = tuple((p, full ^ p) for p in range(8))
Teps = tensor(anti, n)
for maxd in (1, 2):
    if maxd == 1:
        ms = list(matchings(list(range(n)), lambda a, b: ham(a, b) == 1))
    else:
        # matchings with all distances <= 2 and total cost <= 16 are too many; take pairs at distance exactly 2 with p^q in {0011,1100}
        ms = list(matchings(list(range(n)), lambda a, b: (a ^ b) in (3, 12)))
    A = np.array([tensor(m, n) for m in ms], dtype=np.int64).T
    rows = np.nonzero(np.any(A != 0, axis=1) | (Teps != 0))[0]
    M = flint.fmpz_mat(A[rows].tolist())
    aug = flint.fmpz_mat(np.column_stack([A[rows], Teps[rows]]).tolist())
    r1, r2 = M.rank(), aug.rank()
    print(f"D=4 class maxd={maxd}: {len(ms)} matchings, costs {sorted(set(sum(ham(a,b) for a,b in m) for m in ms))}, rank {r1}, with T_eps {r2}, in span: {r1==r2}")
    sys.stdout.flush()
