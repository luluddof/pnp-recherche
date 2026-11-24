# Hill-climb over rank methods Phi = sum_t A_t Fl_{s_t}(.) B_t to maximise R^2/(N^2 c1 c2).
import numpy as np, sys, random
sys.path.insert(0,'.')
from phi import *
N = int(sys.argv[1]); types = [int(c) for c in sys.argv[2]]; m = int(sys.argv[3]); iters = int(sys.argv[4]); seed=int(sys.argv[5])
np.random.seed(seed); random.seed(seed)
D = N*N
def randsparse(r, c, dens):
    M = np.zeros((r,c), dtype=object)
    for i in range(r):
        for j in range(c):
            if random.random() < dens: M[i,j] = random.choice([1,-1,2])
    return M
def score(phi):
    c1,c2,R = profile(phi, N, reps=1)
    if c1 == 0 or c2 == 0: return (R>0)*99, (c1,c2,R)
    return R*R/(D*c1*c2), (c1,c2,R)
phi = [(s, randsparse(m, D, 0.15), randsparse(D, m, 0.15)) for s in types]
best, prof = score(phi)
for it in range(iters):
    t = random.randrange(len(phi)); s,A,B = phi[t]
    A2 = A.copy(); B2 = B.copy()
    M = A2 if random.random()<0.5 else B2
    i = random.randrange(M.shape[0]); j = random.randrange(M.shape[1])
    M[i,j] = 0 if (M[i,j]!=0 and random.random()<0.6) else random.choice([1,-1,2])
    new = phi[:t]+[(s,A2,B2)]+phi[t+1:]
    sc, pr = score(new)
    if sc >= best:
        phi, best, prof = new, sc, pr
print('N',N,'types',types,'m',m,'best ratio %.4f'%best,'(c1,c2,R)=',prof)
