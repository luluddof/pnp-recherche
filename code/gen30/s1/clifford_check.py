# operator-level check (k=2, 4 qubits): Clifford unitaries with symplectic maps A, B (A B = F12) multiply to M12 up to a Pauli
import numpy as np
from sympl_check import blocks
from core import popc
k = 2; n = 2*k; D = 1 << n
def P(a, b):  # Hermitian Pauli i^{a.b} X^a Z^b ; bit j of a <-> qubit j ; vector coords (a bits, b bits)
    z = np.arange(D); M = np.zeros((D, D), complex)
    M[z ^ a, z] = [(-1)**popc(b & int(zz)) for zz in z]
    return (1j)**popc(a & b) * M
def vec2ab(v): a = sum(int(v[j]) << j for j in range(n)); b = sum(int(v[n+j]) << j for j in range(n)); return a, b
def clifford(F, rng):
    gens = []  # images of X_j, Z_j
    for j in range(2*n):
        e = np.zeros(2*n, dtype=np.int64); e[j] = 1; gens.append(P(*vec2ab((F @ e) % 2)))
    M0 = rng.standard_normal((D, D)) + 1j*rng.standard_normal((D, D)); U = np.zeros((D, D), complex)
    for a in range(D):
        for b in range(D):
            img = np.eye(D, dtype=complex)
            for j in range(n):
                if a >> j & 1: img = img @ gens[j]
            for j in range(n):
                if b >> j & 1: img = img @ gens[n+j]
            Ps = np.eye(D, dtype=complex)
            for j in range(n):
                if a >> j & 1: Ps = Ps @ P(1 << j, 0)
            for j in range(n):
                if b >> j & 1: Ps = Ps @ P(0, 1 << j)
            U += img @ M0 @ np.linalg.inv(Ps)
    U /= np.sqrt(np.trace(U.conj().T @ U).real / D)
    return U
rng = np.random.default_rng(0)
F12, FM, Om = blocks(k)
A = np.load('A_k2_3.npy'); B = np.load('B_k2_3.npy')
UA = clifford(A, rng); UB = clifford(B, rng)
print('UA unitary', np.allclose(UA @ UA.conj().T, np.eye(D)), 'UB unitary', np.allclose(UB @ UB.conj().T, np.eye(D)))
# M12 |x,y> = |x+y, y>, x = bits 0..k-1, y = bits k..2k-1
z = np.arange(D); x = z & ((1 << k) - 1); y = z >> k
M12 = np.zeros((D, D)); M12[(x ^ y) | (y << k), z] = 1
W = M12.conj().T @ UA @ UB
# W should be a phase times a Pauli: check its Pauli expansion has a single term
best = max(((a, b) for a in range(D) for b in range(D)), key=lambda ab: abs(np.trace(P(*ab).conj().T @ W)))
c = np.trace(P(*best).conj().T @ W) / D
print('M12^dag UA UB = c * Pauli', best, '|c| =', round(abs(c), 12), 'residual', np.abs(W - c*P(*best)).max())
UM = clifford(FM, rng); z3 = UM @ UM @ UM
print('UM^3 ~ scalar*Pauli? |trace-check|', round(max(abs(np.trace(P(a,b).conj().T @ z3))/D for a in range(D) for b in range(D)), 6))
