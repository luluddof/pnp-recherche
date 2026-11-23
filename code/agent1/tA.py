import random, time
from common import *
import candA_pc as A
rng = random.Random(1)
for name, inst in [('tseitin10', tseitin(10, rng)), ('xor', rand3xor(20, 22, rng)), ('php4', php(4)), ('r30', rand3sat(30, 128, rng))]:
    n, cls = inst
    for d in (3, 4):
        t = time.time(); r = A.refute(n, cls, d); print(name, n, d, truth(n, cls), r, round(time.time()-t, 2))
