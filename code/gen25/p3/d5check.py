# smallest lambda (power of 2 grid) at which (D.5)'s clause u_d >= lambda^{4/5} holds at every active depth
# (paper schedule, child <= 3M/32, M0 = 45.1 n, s = 5.26 n, H = s/floor(sqrt lam)); work with log2 lambda = x
import math
def ok(x, eps):
    lam=2.0**x; C=45.12/5.26
    # active depth d iff C*s*(3/32)^d > s/sqrt(lam)
    dmax=0
    while math.log2(C)+(dmax+1)*math.log2(3/32) > -0.5*x: dmax+=1
    return all(math.log2(eps)+x-d >= 0.8*x for d in range(dmax+1)), dmax
for eps,nm in ((2**-17,'2^-17'),(1/960,'1/960')):
    x=10.0
    while not ok(x,eps)[0]: x+=0.5
    print(f"eps={nm}: (D.5) clause lambda^(4/5)<=u_d holds from log2(lambda) ~ {x} (dmax={ok(x,eps)[1]}), i.e. s >= 2^(2^{x})")
    print(f"   (d=0 alone needs log2 lambda >= {5*math.log2(1/eps):.1f})")
