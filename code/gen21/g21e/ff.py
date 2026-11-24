import numpy as np
def rank_mod(M, p):
    M = np.array(M, dtype=np.int64) % p
    M = M.copy(); r = 0; rows, cols = M.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i, c]: piv = i; break
        if piv is None: continue
        M[[r, piv]] = M[[piv, r]]
        inv = pow(int(M[r, c]), p-2, p)
        M[r] = (M[r]*inv) % p
        nz = np.nonzero(M[:, c])[0]
        for i in nz:
            if i != r: M[i] = (M[i] - M[i, c]*M[r]) % p
        r += 1
        if r == rows: break
    return r
def flat(X, N, S):
    # X indexed [x,y,u,v]; S = tuple of legs for rows
    T = X.reshape(N,N,N,N)
    rest = tuple(i for i in range(4) if i not in S)
    return np.transpose(T, S+rest).reshape(N**len(S), N**len(rest))
