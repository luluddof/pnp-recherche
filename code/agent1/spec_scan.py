"""Smallest clause count m at which spectral refutation (D) certifies UNSAT, vs n."""
import random
from common import rand3sat
from candBCD import spectral_refute
rng = random.Random(5)
for n in (30, 60, 100, 150):
    m = 4 * n
    while True:
        oks = sum(spectral_refute(*rand3sat(n, m, rng))[0] for _ in range(3))
        if oks >= 2: break
        m = int(m * 1.25)
    print(f'n={n} m*={m} m*/n={m/n:.1f} m*/n^1.5={m/n**1.5:.2f} (bound/m at 4.26n = {spectral_refute(*rand3sat(n, round(4.26*n), rng))[1]:.2f})', flush=True)
