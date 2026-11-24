import numpy as np
from scipy.optimize import least_squares
from distlm import res, unpack
from dist import swap
N,p,q=3,2,2; T=swap(N).astype(complex); rng=np.random.default_rng(11)
for r in range(3):
    s = least_squares(res, rng.standard_normal(4*(p+q)*N*N)*0.7, args=(T,N,p,q), method='trf', max_nfev=3000, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    a,b,G,H = unpack(s.x,N,p,q)
    A = np.einsum('kxy,kuv->xyuv', a, b); B = np.einsum('kxu,kyv->xyuv', G, H); E = T - A - B
    sv = lambda X, perm: np.round(np.linalg.svd(X.transpose(perm).reshape(9,9), compute_uv=False), 4)
    print('dist^2', round(2*s.cost, 6))
    print('  residual  Fl1 sv', sv(E,(0,1,2,3)), ' Fl2 sv', sv(E,(0,2,1,3)), ' Fl3 sv', sv(E,(0,3,1,2)))
    print('  A Fl1 sv', sv(A,(0,1,2,3))[:3], ' B Fl2 sv', sv(B,(0,2,1,3))[:3], ' |A|,|B|', round(np.linalg.norm(A),3), round(np.linalg.norm(B),3))
    print('  B Fl1 sv', sv(B,(0,1,2,3)))
    print('  A+B Fl1 sv', sv(A+B,(0,1,2,3)))
