# G24e: Levenberg-Marquardt on the FULL polynomial system  X = U V^T + sum_l G_l (x) H_l   (complex, holomorphic residual)
# r[(x,y),(u,v)] = X - sum_k U[(x,y),k] V[(u,v),k] - sum_l G[l,x,u] H[l,y,v].
# Optional Tikhonov term mu*||z||^2 (continuation in mu: exact point => residual -> 0 with bounded norm; border => norm blows up).
import numpy as np, sys, json, time

def swapmat(N):
    X = np.zeros((N, N, N, N))
    for x in range(N):
        for y in range(N): X[x, y, y, x] = 1
    return X.reshape(N * N, N * N)

def unpack(z, N, p, q):
    s = N * N; i = 0
    U = z[i:i + s * p].reshape(s, p); i += s * p
    V = z[i:i + s * p].reshape(s, p); i += s * p
    G = z[i:i + q * s].reshape(q, N, N); i += q * s
    H = z[i:i + q * s].reshape(q, N, N)
    return U, V, G, H

def resid(z, N, p, q, X):
    U, V, G, H = unpack(z, N, p, q)
    C = np.einsum('lxu,lyv->xyuv', G, H).reshape(N * N, N * N)
    return X - U @ V.T - C

def jac(z, N, p, q):
    U, V, G, H = unpack(z, N, p, q); s = N * N; I = np.eye(s); In = np.eye(N)
    JU = -np.einsum('ac,bk->abck', I, V).reshape(s * s, s * p)
    JV = -np.einsum('ak,bc->abck', U, I).reshape(s * s, s * p)
    JG = -np.einsum('xa,ub,lyv->xyuvlab', In, In, H).reshape(s * s, q * s)
    JH = -np.einsum('lxu,ya,vb->xyuvlab', G, In, In).reshape(s * s, q * s)
    return np.hstack([JU, JV, JG, JH])

def lm(z, N, p, q, mu=0.0, iters=200, lam=1e-3, verbose=False, tol=1e-28):
    X = swapmat(N).astype(z.dtype)
    hist = []
    def cost(z):
        r = resid(z, N, p, q, X); return float(np.vdot(r, r).real), float(np.vdot(r, r).real + mu * np.vdot(z, z).real)
    f0, F0 = cost(z)
    for it in range(iters):
        r = resid(z, N, p, q, X).ravel(); J = jac(z, N, p, q)
        A = J.conj().T @ J; g = J.conj().T @ r
        if mu: A = A + mu * np.eye(len(z)); g = g + mu * z
        dA = np.diag(A).real
        while True:
            try: dz = np.linalg.solve(A + lam * np.diag(dA + 1e-12), -g)
            except np.linalg.LinAlgError: lam *= 10; continue
            f1, F1 = cost(z + dz)
            if F1 < F0: z = z + dz; f0, F0 = f1, F1; lam = max(lam / 3, 1e-15); break
            lam *= 4
            if lam > 1e16: break
        hist.append((it, f0, float(np.abs(z).max()), float(np.linalg.norm(z))))
        if verbose and (it % 10 == 0 or f0 < tol): print(f"  it{it} d2={f0:.3e} max|z|={hist[-1][2]:.3f} |z|={hist[-1][3]:.2f} lam={lam:.1e}", flush=True)
        if f0 < tol or lam > 1e16: break
    return z, hist

def from_vp(zvp, N, p, q):
    """vp.py point (B side has q terms, real/imag stacked) -> full LM variables (A = best rank-p approx of F - B)."""
    n = len(zvp) // 2; w = zvp[:n] + 1j * zvp[n:]; s = N * N
    G = w[:q * s].reshape(q, N, N); H = w[q * s:].reshape(q, N, N)
    M = swapmat(N) - np.einsum('lxu,lyv->xyuv', G, H).reshape(s, s)
    Us, S, Vh = np.linalg.svd(M)
    U = Us[:, :p] * np.sqrt(S[:p]); V = (Vh[:p, :].T) * np.sqrt(S[:p])
    return np.concatenate([U.ravel(), V.ravel(), G.ravel(), H.ravel()])
