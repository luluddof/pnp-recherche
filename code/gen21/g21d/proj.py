# dims of coordinate projections of J(2,2) (N=3) to sub-formats (restrict each leg to a subset of indices)  [E, Terracini mod PR]
import numpy as np, itertools
from common import *
from terracini import jac
rng = np.random.default_rng(2); Jc = jac(3, 2, 2, rng).reshape(-1, 3, 3, 3, 3)
for fmt in [(2,3,3,3),(3,2,3,3),(3,3,2,3),(3,3,3,2),(2,3,3,2),(2,2,3,3),(2,3,2,3),(3,2,2,3),(2,2,2,2),(2,3,3,1),(1,3,3,3)]:
    sub = Jc[:, :fmt[0], :fmt[1], :fmt[2], :fmt[3]].reshape(len(Jc), -1)
    print(fmt, 'ambient', sub.shape[1], 'dim proj J(2,2) =', rank(sub), flush=True)
