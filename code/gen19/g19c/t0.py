import numpy as np, time, sys
from starlib import *
inst = Inst(2,2,2)
print(inst.d, inst.K, inst.D, inst.npar, inst.prod_r)
# check P vs hexlib
sys.path.insert(0,'/home/user/pnp-recherche/code/gen18/g18f')
import hexlib
print('P matches hexlib:', np.array_equal(inst.P, hexlib.perm_P(2)))
rng = np.random.default_rng(0)
x = rng.normal(size=inst.npar); T = inst.unpack(x)
# to3/tomat inverse; jvec/jtvec adjoint; jvec = derivative
M = rng.normal(size=(inst.D,inst.D)); print('to3/tomat', np.abs(inst.tomat(inst.to3(M))-M).max())
dx = rng.normal(size=inst.npar); G = rng.normal(size=inst.K)
print('adjoint', np.sum(inst.jvec3(T, inst.unpack(dx))*G), inst.jtvec3(T,G)@dx)
e=1e-6; print('deriv', np.abs((inst.build3(inst.unpack(x+e*dx))-inst.build3(inst.unpack(x-e*dx)))/(2*e)-inst.jvec3(T,inst.unpack(dx))).max())
# compare build with hexlib build for star0
T0 = T[0:3]; d=4
B0 = T0[0].reshape(d,d,4,4); B1 = T0[1].reshape(d,d,4,1); B2 = T0[2].reshape(d,d,1,4)
Mh = hexlib.build(B0,B1,B2)  # 'aiAG,bjAB,ckBG': B0[.,.,al,ga], B1[.,.,al,be], B2[.,.,be,ga]
print('star0 vs hexlib', np.abs(inst.tomat(inst.star3(T,0))-Mh).max())
print('regauge', np.abs(inst.build3(inst.unpack(inst.regauge(x)))-inst.build3(T)).max())
for r in [3,2,1,0]:
    fs=[]; t0=time.time()
    for seed in range(10):
        x0 = inst.random_x(np.random.default_rng(seed), 0.1)
        x, f, h = lm_solve(inst, r, x0, maxit=300)
        fs.append(f)
    nm = norms(inst, x, r)
    print(f"N=2 r={r}: f={np.sort(fs)}  its={len(h)} time={time.time()-t0:.1f}s  last norms={nm}")
