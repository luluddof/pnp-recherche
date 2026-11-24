# r-colored partitions: L(E) = r*phi(E/r) asymptotically; ratio phi(2x)/phi(x) vs sqrt2
import numpy as np
from scipy.optimize import minimize_scalar
def F(t):
    k=np.arange(1,200000)
    return -np.sum(np.log1p(-np.exp(-k*t)))
def phi(x):
    r=minimize_scalar(lambda lt: F(np.exp(lt))+x*np.exp(lt),bounds=(-12,6),method='bounded')
    return r.fun
for x in [1e-4,1e-3,1e-2,0.03,0.1,0.3,1,3,10,30,100,1000]:
    print(x, phi(2*x)/phi(x))
