# [N] sanity check of every inequality in the proof of Theorem A (c(N;p,1) <= p + |T| + gamma^2), random + adversarial inputs.
import numpy as np
rng = np.random.default_rng(7)
def swap(N):
    F = np.zeros((N*N, N*N))
    for x in range(N):
        for y in range(N): F[x*N+y, y*N+x] = 1
    return F
def realign(X, N):
    return X.reshape(N, N, N, N).transpose(0, 2, 1, 3).reshape(N*N, N*N)
def cplx(*s): return rng.standard_normal(s) + 1j*rng.standard_normal(s)
worst = 0; viol = 0
for trial in range(3000):
    N = rng.integers(2, 5); p = rng.integers(1, 4); F = swap(N)
    G, Gp = cplx(N, N), cplx(N, N)
    if trial % 3 == 0:   # near-extremal: low-rank G, G' and A antisymmetric-ish
        W = np.linalg.qr(cplx(N, 2))[0]; G = W @ W.conj().T + 0.1*cplx(N, N); Gp = G.T + 0.1*cplx(N, N)
    a, b = cplx(N*N, p), cplx(N*N, p)
    if trial % 3 == 0:
        a = (np.eye(N*N) - F) @ a; b = a.copy()
    A = a @ b.conj().T * rng.choice([0.3, 1, 3])
    B = np.kron(G, Gp); X = A + B
    # P = projection onto row space of A
    U_, s_, Vh = np.linalg.svd(A); r = (s_ > 1e-10*s_[0]).sum(); P = Vh[:r].conj().T @ Vh[:r]; Q = np.eye(N*N) - P
    # product singular terms
    u1, s1, v1 = np.linalg.svd(G); u2, s2, v2 = np.linalg.svd(Gp)
    terms = [(s1[i]*s2[k], np.kron(u1[:, i], u2[:, k]), np.kron(v1[i].conj(), v2[k].conj())) for i in range(N) for k in range(N)]
    th = 0.5 + 1e-9
    T = [j for j, (be, w, z) in enumerate(terms) if np.linalg.norm(P @ z)**2 >= th]
    BT = sum((terms[j][0]*np.outer(terms[j][1], terms[j][2].conj()) for j in T), np.zeros((N*N, N*N), complex))
    BR = B - BT
    # checks
    c1 = abs(np.trace(F @ BT @ Q)) <= np.sqrt(len(T))*np.linalg.norm(BT @ Q) + 1e-9
    c2 = np.linalg.norm(BR @ P)**2 <= th*np.linalg.norm(BR)**2 + 1e-9
    c3 = abs(np.trace(F @ BR)) <= np.linalg.svd(realign(BR, N), compute_uv=False).sum() + 1e-9
    rk = np.linalg.matrix_rank(realign(BR, N), tol=1e-8)
    c4 = rk <= 1 + len(T)
    c5 = abs(np.linalg.norm(B @ Q)**2 - np.linalg.norm(BT @ Q)**2 - np.linalg.norm(BR @ Q)**2) < 1e-8*max(1, np.linalg.norm(B)**2)
    gam2 = (np.sqrt(1 + len(T)) + np.sqrt(p*th))**2/(1 - th)
    ratio = abs(np.trace(F @ X))**2/np.linalg.norm(X)**2
    c6 = ratio <= p + len(T) + gam2 + 1e-9
    if not (c1 and c2 and c3 and c4 and c5 and c6): viol += 1; print("VIOLATION", trial, c1, c2, c3, c4, c5, c6)
    worst = max(worst, ratio/(2*p + max(m for m in range(1, 20) if m*(m-1)//2 <= p)))
print("violations:", viol, " max ratio/(2p+m) over samples:", worst)
