import random, time
from common import *
from candBCD import *
rng = random.Random(2)
n, cls = rand3sat(50, 213, rng); print('truth', truth(n, cls))
t=time.time(); print('clusterLP feasible', cluster_lp_feasible(n, cls), time.time()-t)
t=time.time(); a=candB(n, cls, rng); print('B', a is not None and check(cls, a), time.time()-t)
t=time.time(); a=bp_decimate(n, cls, rng); print('C', a is not None and check(cls, a), time.time()-t)
t=time.time(); print('D', spectral_refute(n, cls), time.time()-t)
for m in [2000, 8000, 20000, 60000]:
    n, cls = rand3sat(40, m, rng); print('D n=40 m=',m, spectral_refute(n, cls))
