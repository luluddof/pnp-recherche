import numpy as np, random
from core import *
from verify_iso import XZ, sigma
random.seed(2)
for k in (1, 2):
    N = 2**k; S = sigma(k); Mm = Mperm(k)
    paulis = [(a1,a2,c,d) for a1 in range(N) for a2 in range(N) for c in range(N) for d in range(N)]
    coef = {}
    for s in paulis:
        P = pauli_M(k, *s); coef[s] = int(np.sum(P * Mm))  # Tr(P^T M) ; P real orthogonal-ish: P^dagger = P^T
    assert all(abs(v) == 1 for v in coef.values()), set(coef.values())
    for trial in range(8):
        sub = random.sample(paulis, random.randint(N**4//2, N**4))
        C = sum(coef[s]*np.kron(np.kron(XZ(k,s[0],s[2]), XZ(k,s[0]^s[1],s[2]^s[3])), XZ(k,s[1],s[3])) for s in sub)
        Cp = sum(coef[s]*pauli_M(k,*s) for s in sub)
        r3 = rank2(N*N*S - C); r2 = rank2(N*N*Mm - Cp)
        print(k, len(sub), r3, r2, 'OK' if r3[0] == N*r2[0] else 'MISMATCH')
