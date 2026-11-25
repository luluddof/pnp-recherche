# Claim 1 + 2 checks. Real Paulis P_ab = X^a Z^b on k qubits; K = sum kappa_ab P(x)P.
import numpy as np, itertools, sys
from math import comb
def pauli(a, b, k):
    N = 2**k; M = np.zeros((N, N))
    for y in range(N):
        sgn = (-1)**bin(b & y).count('1')   # Z^b first, then X^a
        M[y ^ a, y] = sgn
    return M
def swap(N):
    n = N*N; S = np.zeros((n, n)); x = np.repeat(np.arange(N), N); y = np.tile(np.arange(N), N); S[y*N+x, x*N+y] = 1; return S
def bell_count(k, kap):
    # phi(c,d) = sum kappa_ab (-1)^{a.d + b.c}; beta = (-1)^{c.d}
    N = 2**k; cnt = 0
    for c in range(N):
        for d in range(N):
            phi = sum(v * (-1)**(bin(a & d).count('1') + bin(b & c).count('1')) for (a, b), v in kap.items())
            if abs(phi - (-1)**bin(c & d).count('1')) > 1e-9: cnt += 1
    return cnt
rng = np.random.default_rng(7)
# claim 1: plain truncation, random signs/phases of modulus 1/N, q < N -> full rank (qudit Weyl too)
for N in [3, 4, 5, 8]:
    w = np.exp(2j*np.pi/N); Xs = np.roll(np.eye(N), 1, 0); Zs = np.diag(w**np.arange(N))
    W = [np.linalg.matrix_power(Xs, a) @ np.linalg.matrix_power(Zs, b) for a in range(N) for b in range(N)]
    full = sum(np.kron(P, P.conj().T) for P in W) / N
    assert np.allclose(full, swap(N))
    worst = N*N
    for t in range(30):
        S = rng.choice(N*N, N-1, replace=False); ph = np.exp(2j*np.pi*rng.random(N-1))
        K = sum(ph[i]*np.kron(W[s], W[s].conj().T) for i, s in enumerate(S)) / N
        worst = min(worst, np.linalg.matrix_rank(swap(N) - K))
    print(f"claim1 N={N}: Weyl expansion of SWAP ok; min rank(SWAP-K) over 30 random S (q=N-1) = {worst} (N^2={N*N})")
# claim 2: Bell formula == explicit rank
for k in [2, 3]:
    N = 2**k; ok = 0
    for t in range(40):
        q = rng.integers(1, N*N); S = rng.choice(N*N, q, replace=False)
        kap = {(s // N, s % N): float(rng.integers(-3, 4)) / rng.integers(1, 4) for s in S}
        if t % 2: kap = {key: (1.0 if v >= 0 else -1.0) / N for key, v in kap.items()}
        K = sum(v * np.kron(pauli(a, b, k), pauli(a, b, k)) for (a, b), v in kap.items())
        ok += (np.linalg.matrix_rank(swap(N) - K) == bell_count(k, kap))
    print(f"claim2 N={N}: Bell count == explicit rank in {ok}/40 random cases")
