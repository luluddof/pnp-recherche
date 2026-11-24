# Reynolds projection onto K-invariants, K a finite subgroup of Stab(F) inside the symmetry group of (sigma_p, tau_q):
#   X -> (a(x)b) X (b^-1 (x) a^-1),  X -> X^T (Fl1 transpose),  X -> F X F.   pi(J) = join(pi(sigma_p), pi(tau_q)) (pi linear);
#   F notin closure pi(J)  ==>  F notin J.  Here: dims (Terracini, mod PR) of W^K and pi(J).
import numpy as np, itertools, sys
from common import *
from terracini import jac
N = 3
def op_ab(a, b):
    ai = np.round(np.linalg.inv(a)).astype(np.int64); bi = np.round(np.linalg.inv(b)).astype(np.int64)
    M = np.zeros((81,81), dtype=np.int64)
    for i in range(81):
        X = np.zeros(81, dtype=np.int64); X[i] = 1; X = X.reshape(9,9)
        Y = np.kron(a,b) @ X @ np.kron(bi, ai); M[:, i] = Y.ravel()
    return M
def op_T():
    M = np.zeros((81,81), dtype=np.int64)
    for i in range(81):
        X = np.zeros(81, dtype=np.int64); X[i] = 1; M[:, i] = X.reshape(9,9).T.ravel()
    return M
Fm = swap(N).reshape(9,9)
def op_F():
    M = np.zeros((81,81), dtype=np.int64)
    for i in range(81):
        X = np.zeros(81, dtype=np.int64); X[i] = 1; M[:, i] = (Fm @ X.reshape(9,9) @ Fm).ravel()
    return M
def closure(gens):
    G = {np.eye(81, dtype=np.int64).tobytes(): np.eye(81, dtype=np.int64)}; frontier = list(G.values())
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                k = h @ g; key = k.tobytes()
                if key not in G: G[key] = k; new.append(k)
        frontier = new
        if len(G) > 5000: raise Exception('too big')
    return list(G.values())
def perm(s):
    P = np.zeros((3,3), dtype=np.int64)
    for i, j in enumerate(s): P[j, i] = 1
    return P
I3 = np.eye(3, dtype=np.int64); cyc = perm([1,2,0]); tr = perm([1,0,2]); sg = np.diag([-1,1,1]).astype(np.int64)
def analyse(name, gens, p=2, q=2):
    K = closure(gens); R = sum(K) % PR   # |K| * Reynolds
    dW = rank(R); rng = np.random.default_rng(1)
    Jc = jac(N, p, q, rng)              # rows = tangent vectors of the parametrisation
    dJ = rank((Jc @ R.T) % PR)
    Fin = rank(np.vstack([(Jc @ R.T) % PR, Fm.ravel()[None, :]]))  # is F in the tangent span? (only a hint)
    print(f"{name}: |K|={len(K)} dim W^K={dW} dim pi(J({p},{q}))>={dJ}  {'PROPER' if dJ < dW else 'fills'}", flush=True)
    return K, R
if __name__ == "__main__":
    tests = {
      'S3 diag (a=b)': [op_ab(cyc,cyc), op_ab(tr,tr)],
      'S3xS3': [op_ab(cyc,I3), op_ab(tr,I3), op_ab(I3,cyc), op_ab(I3,tr)],
      'S3 diag + T + AdF': [op_ab(cyc,cyc), op_ab(tr,tr), op_T(), op_F()],
      'Z3 diag': [op_ab(cyc,cyc)],
      'Z3 diag + T + AdF': [op_ab(cyc,cyc), op_T(), op_F()],
      'T + AdF': [op_T(), op_F()],
      'signs a (x) 1': [op_ab(sg,I3), op_ab(cyc@sg@cyc.T,I3)],
      'S3 on a only': [op_ab(cyc,I3), op_ab(tr,I3)],
      'S3 on a only + T + AdF': [op_ab(cyc,I3), op_ab(tr,I3), op_T(), op_F()],
      'S3xS3 + T + AdF': [op_ab(cyc,I3), op_ab(tr,I3), op_ab(I3,cyc), op_ab(I3,tr), op_T(), op_F()],
      'signed S3 diag': [op_ab(cyc,cyc), op_ab(tr,tr), op_ab(sg,sg)],
      'signed S3 diag + T + AdF': [op_ab(cyc,cyc), op_ab(tr,tr), op_ab(sg,sg), op_T(), op_F()],
    }
    for k, g in tests.items(): analyse(k, g)
