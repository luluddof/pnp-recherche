# Banded factors (bandwidth w): K supported on ((i,k),(j,l)) with |i-j|<=w, |k-l|<=w. Explicit K achieving (N-w)(N-w-1)/2.
import numpy as np
from modrank import rank_mod
for N in [5, 6, 8, 10]:
    for w in range(0, N):
        NN = N*N; S = np.zeros((NN, NN), dtype=np.int64); K = np.zeros_like(S)
        for i in range(N):
            for k in range(N):
                S[k*N+i, i*N+k] = 1
                if abs(i-k) <= w: K[k*N+i, i*N+k] = 1          # free: cancel
                elif i < k:                                     # 2x2 block [[ -1, 1],[1,-1]] -> rank 1 (K diag = 1)
                    K[i*N+k, i*N+k] = 1; K[k*N+i, k*N+i] = 1
        # check support
        ok = all(abs(a//N - b//N) <= w and abs(a % N - b % N) <= w for a, b in zip(*np.nonzero(K)))
        r = rank_mod(S - K, 1000003)
        print(N, w, "support ok", ok, "rank", r, "improved bound", (N-w)*(N-w-1)//2, "C2 bound", max(0,(N-2*w)*(N-2*w-1)//2))
