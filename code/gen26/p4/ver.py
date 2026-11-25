# independent exact check: signed-partial-permutation property + g19b/verify.py (flint rank over Q, explicit factorisation)
import sys, json, numpy as np
import verify_g19b as verify
for fn in sys.argv[1:]:
    d = json.load(open(fn)); N = d['N']
    ok = all(set(np.unique(np.array(m))) <= {0, 1, -1} and (np.abs(np.array(m)).sum(0) <= 1).all() and (np.abs(np.array(m)).sum(1) <= 1).all() for m in d['E'] + d['C'])
    print(fn, "signed partial perms:", ok, " q =", len(d['E']))
    verify.verify(N, d['E'], d['C'], fn)
