# [N] kappa(D) = sup_H tr(Y conj(M)^{-1} Y^*),  Y = H - D H^T D,  M = ||H||^2 I - D H H^* D,  D = diag(sigma), sum sigma^2 = 1.
# (= kappa(P) for P = |b><b|, b with Schmidt coefficients sigma; see notes section 1.)
import numpy as np, sys
from scipy.optimize import minimize

def val(H, D):
    H = H/np.linalg.norm(H)
    Y = H - D @ H.T @ D
    M = np.eye(len(D)) - D @ H @ H.conj().T @ D
    return float(np.real(np.trace(Y @ np.linalg.solve(M.conj(), Y.conj().T))))

def kapD(sig, R=8, seed=0):
    N = len(sig); D = np.diag(sig); rng = np.random.default_rng(seed); best = -1; bh = None
    for r in range(R):
        f = lambda z: -val((z[:N*N] + 1j*z[N*N:]).reshape(N, N), D)
        s = minimize(f, rng.standard_normal(2*N*N), method='BFGS', options=dict(gtol=1e-11, maxiter=4000))
        if -s.fun > best: best, bh = -s.fun, s.x
    H = (bh[:N*N] + 1j*bh[N*N:]).reshape(N, N)
    return best, H/np.linalg.norm(H)

if __name__ == "__main__":
    np.set_printoptions(precision=3, suppress=True, linewidth=150)
    for sig in [(1,0), (np.cos(.3), np.sin(.3)), (np.cos(.6), np.sin(.6)), (1,1), (1,1,0), (1,1,1), (2,1,1), (3,2,1), (1,1,1,1)]:
        s = np.array(sig, float); s /= np.linalg.norm(s)
        k, H = kapD(s)
        print(np.round(s,4), "kappa=%.8f" % k, " 2+2s1s2=%.6f" % (2 + 2*s[0]*s[1]))
