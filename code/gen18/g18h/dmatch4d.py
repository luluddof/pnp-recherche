"""D=4, lambda=1: ALL 10064 eps-matchings of total cost <= 10, restricted to 3000 random balanced coordinates, mod p.
'not in span' on a coordinate subset mod p  =>  not in span mod p (restriction is linear).  [Np]"""
import numpy as np, flint, sys, time, random
from dmatch import tensor
from dmatch4b_gen import gen
p = 2**31 - 1; n = 16
anti = tuple((q, 15 ^ q) for q in range(8)); Teps = tensor(anti, n)
bal = [i for i in range(1 << 16) if bin(i).count("1") == 8]
random.seed(7); sub = sorted(random.sample(bal, 3000))
t0 = time.time()
ms = list(gen(list(range(16)), 10))
A = np.array([tensor(m, n)[sub] for m in ms], dtype=np.int64) % p
r1 = flint.nmod_mat(A.tolist(), p).rank()
r2 = flint.nmod_mat(np.vstack([A, Teps[sub] % p]).tolist(), p).rank()
print(f"C=10: all {len(ms)} matchings on 3000 coordinates: rank mod p {r1}, with T_eps {r2}, in span: {r1 == r2} ({time.time()-t0:.0f}s)")
