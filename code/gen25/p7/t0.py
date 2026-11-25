import numpy as np
from spread import *; from load import *
R, C, p, q = load_g19('exact_5_6_6.json'); N = R.shape[0]
d = analyze(R, C); print(fmt(d), p, q)
rng = np.random.default_rng(1); a = rng.standard_normal((N,N)); b = rng.standard_normal((N,N))
R2, C2 = act(R, a, b), act(C, a, b); print('F fixed', np.abs(act(swapT(N), a, b) - swapT(N)).max())
print(fmt(analyze(R2, C2)))
# gradient check
m, beta = 24, 30
f, G4R, G4C, Lt = objective(R, C, p, q, m, beta)
gXR, gYR = grads_XY(G4R, R); gXC, gYC = grads_XY(G4C, C)
X = rng.standard_normal((N,N)); Y = rng.standard_normal((N,N)); e = 1e-6
from scipy.linalg import expm
fp = objective(act(R, expm(e*X), expm(e*Y)), act(C, expm(e*X), expm(e*Y)), p, q, m, beta)[0]
print('fd', (fp - f)/e, 'an', ((gXR+gXC)*X).sum() + ((gYR+gYC)*Y).sum())
