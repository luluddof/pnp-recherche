# check rank(sigma - C) = N rank(M - C') for random sparse Pauli combos, k=1,2 (integer coefficients)
import numpy as np, itertools, random
from core import *
def kron(*m):
    r = np.array([[1]], dtype=np.int64)
    for x in m: r = np.kron(r, x)
    return r
X1 = np.array([[0,1],[1,0]]); Z1 = np.array([[1,0],[0,-1]]); I1 = np.eye(2, dtype=np.int64)
def XZ(k, a, c):  # X^a Z^c on k qubits, bit j <-> qubit j; basis index = integer with bit j = qubit j
    n = 2**k; P = np.zeros((n, n), dtype=np.int64)
    for i in range(n): P[i ^ a, i] = (-1)**pc(c & i)
    return P
def sigma(k):
    N = 2**k; S = np.zeros((N**3, N**3), dtype=np.int64)
    for u1 in range(N):
        for u2 in range(N):
            for u3 in range(N):
                S[(u3*N + u1)*N + u2, (u1*N + u2)*N + u3] = 1
    return S
random.seed(1)
for k in (1, 2):
    N = 2**k; S = sigma(k); Mm = Mperm(k)
    for trial in range(6):
        m = random.randint(1, 3*N*N)
        C = np.zeros_like(S); Cp = np.zeros_like(Mm)
        for _ in range(m):
            a1, a2, c, d = [random.randrange(N) for _ in range(4)]; w = random.randint(-3, 3)
            C += w * np.kron(np.kron(XZ(k, a1, c), XZ(k, a1 ^ a2, c ^ d)), XZ(k, a2, d))
            Cp += w * pauli_M(k, a1, a2, c, d)
        for t in (1, 2):  # scale so that sigma - C/t vs M - C'/t
            r3 = rank2(t*S - C); r2 = rank2(t*Mm - Cp)
            print(k, trial, m, t, r3, r2, 'OK' if r3[0] == N*r2[0] else 'MISMATCH')
