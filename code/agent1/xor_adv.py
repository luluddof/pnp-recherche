"""Adversarial: planted (hence SAT) random 3-XOR near threshold, for SAT-side heuristics B, C,
versus Gaussian elimination (exact, poly)."""
import random, sys
from common import *
from candBCD import *
rng = random.Random(99)
for n in (100, 200, 400):
    okB = okC = 0; T = 6
    for _ in range(T):
        plant = [None] + [rng.random() < .5 for _ in range(n)]
        cls = []
        for _ in range(round(0.85 * n)):
            vs = rng.sample(range(1, n + 1), 3)
            cls += xor_clauses(vs, sum(plant[v] for v in vs) % 2)
        assert check(cls, plant)
        a = candB(n, cls, rng); okB += a is not None and check(cls, a)
        a = bp_decimate(n, cls, rng); okC += a is not None and check(cls, a)
    print('planted 3-XOR density .85 n', n, 'B(LP+WalkSAT 20n^2)', f'{okB}/{T}', 'C(BP)', f'{okC}/{T}', flush=True)
