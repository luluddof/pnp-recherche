"""Adversarial scaling of the SAT-side heuristics B (LP+WalkSAT, 20n^2 flips) and C (BP decimation)
at density 4.26 for large n. Ground truth: CaDiCaL with a conflict budget (unknowns dropped)."""
import random, time, sys
from common import *
from candBCD import *
from pysat.solvers import Cadical153
rng = random.Random(777)
dens = float(sys.argv[1]) if len(sys.argv) > 1 else 4.26
for n in (200, 400, 800):
    cnt = [0, 0, 0, 0]; tB = tC = 0
    for _ in range(10):
        n_, cls = rand3sat(n, round(dens * n), rng)
        with Cadical153(bootstrap_with=cls) as s:
            s.conf_budget(3000000); r = s.solve_limited()
        if r is None: cnt[3] += 1; continue
        if not r: cnt[2] += 1; continue
        cnt[0] += 1
        t = time.time(); a = candB(n, cls, rng); tB += time.time() - t; b = a is not None and check(cls, a)
        t = time.time(); a = bp_decimate(n, cls, rng); tC += time.time() - t; c = a is not None and check(cls, a)
        cnt[1] += 0; print(n, 'SAT', 'B', b, 'C', c, flush=True)
    print('n', n, 'SAT', cnt[0], 'UNSAT', cnt[2], 'unknown', cnt[3], 'tB', round(tB, 1), 'tC', round(tC, 1), flush=True)
