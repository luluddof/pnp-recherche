import sys; sys.path.insert(0,'.')
from validate import cd_terms
from ev import *
for n in range(1,6):
    N=1<<n; print(N, best_scale(N, cd_terms(n)))
