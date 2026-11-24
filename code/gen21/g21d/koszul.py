# Flattenings and Koszul (Young) flattenings of 4-leg tensors: generic rank on J(p,q) (mod PR, a lower bound for the
# char-0 generic rank) versus EXACT rank at the swap (over Z via flint fmpz_mat).  Separation needs rank(F) > generic rank on J.
import numpy as np, itertools, sys
from flint import fmpz_mat
from common import *
def kosz(X, Agrp, S, k, N):
    legs = [0,1,2,3]; Sc = [l for l in legs if l not in Agrp and l not in S]
    Y = X.transpose(list(Agrp)+list(S)+list(Sc))
    nA = N**len(Agrp); nS = N**len(S); nSc = N**len(Sc)
    Y = Y.reshape(nA, nS, nSc)
    Ik = list(itertools.combinations(range(nA), k)); Ik1 = list(itertools.combinations(range(nA), k+1))
    pos = {J:i for i,J in enumerate(Ik1)}
    M = np.zeros((len(Ik1)*nSc, len(Ik)*nS), dtype=np.int64)
    for ci, I in enumerate(Ik):
        for a in range(nA):
            if a in I: continue
            J = tuple(sorted(I + (a,))); sg = (-1)**sum(1 for i in I if i < a); ri = pos[J]
            # M[(ri, s'), (ci, s)] += sg * Y[a, s, s']
            M[ri*nSc:(ri+1)*nSc, ci*nS:(ci+1)*nS] += sg * Y[a].T
    return M
def zrank(M):
    return fmpz_mat(M.shape[0], M.shape[1], [int(v) for v in M.ravel()]).rank()
if __name__ == "__main__":
    N, p, q = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rng = np.random.default_rng(5); F = swap(N)
    J = [ptJ(N,p,q,rng) for _ in range(2)]
    s1 = [ptJ(N,1,0,rng) for _ in range(2)]; s2 = [ptJ(N,0,1,rng) for _ in range(2)]
    res = []
    groups = [(a,) for a in range(4)] + [g for g in itertools.combinations(range(4), 2)]
    for Ag in groups:
        rest = [l for l in range(4) if l not in Ag]
        nA = N**len(Ag)
        for k in range(0, min(nA, 3 if len(Ag)==2 else nA)):
            for r in range(0, len(rest)+1):
                for S in itertools.combinations(rest, r):
                    if k == 0 and r == 0 and len(Ag)==1: pass
                    shape = kosz(F, Ag, S, k, N).shape
                    if shape[0]*shape[1] > 1300*1300: continue
                    rF = zrank(kosz(F, Ag, S, k, N))
                    rJ = max(rank(kosz(X, Ag, S, k, N)) for X in J)
                    c1 = max(rank(kosz(X, Ag, S, k, N)) for X in s1); c2 = max(rank(kosz(X, Ag, S, k, N)) for X in s2)
                    tag = "SEPARATES" if rF > rJ else ""
                    line = f"A={Ag} k={k} S={S} shape={shape} c1={c1} c2={c2} rank_J>={rJ} rank_F={rF} {tag}"
                    print(line, flush=True)
