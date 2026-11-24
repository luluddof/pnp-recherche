"""D=4, lambda=1, eps-matching class, mod p = 2^31-1 (NOT exact over Q: [N mod p]).
Smallest C with T_eps in span{eps-matchings of total cost <= C} (cost = sum of Hamming distances; floor 8, full 32)."""
import numpy as np, flint, sys, time, random
from dmatch import tensor
from dmatch4b_gen import gen
p = 2**31 - 1
n = 16; full = 15
anti = tuple((q, full ^ q) for q in range(8))
Teps = tensor(anti, n)
bal = [i for i in range(1 << 16) if bin(i).count("1") == 8]
random.seed(1)
for C in (10, 12, 14):
    t0 = time.time()
    ms = list(gen(list(range(16)), C)) if C <= 12 else None
    if ms is None: break
    tot = len(ms)
    if len(ms) > 4000: ms = random.sample(ms, 4000)
    A = np.array([tensor(m, n)[bal] for m in ms], dtype=np.int64) % p
    M = flint.nmod_mat(A.tolist(), p); r1 = M.rank()
    aug = flint.nmod_mat(np.vstack([A, Teps[bal] % p]).tolist(), p); r2 = aug.rank()
    print(f"C={C}: {tot} matchings, used {len(ms)}; rank mod p {r1}, with T_eps {r2}, in span: {r1 == r2} ({time.time()-t0:.0f}s)"); sys.stdout.flush()
