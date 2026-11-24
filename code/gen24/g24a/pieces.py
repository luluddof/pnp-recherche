# per-piece (rr, rc) profiles of Hrubes's C-indexed pieces and threshold trade-offs, exact
import numpy as np, scipy.sparse as sp, sys, time, json
from math import comb
import ev, hr
def rr(T): N = T.shape[0]; return ev.rank(sp.csr_matrix(T.reshape(N*N, N*N)))[0]
def rc(T): N = T.shape[0]; return ev.rank(sp.csr_matrix(T.transpose(0, 2, 1, 3).reshape(N*N, N*N)))[0]
def piece_tensors(S, ci):
    """returns (P1, P2) as 4-tensors T[x=A, y=i, u=j, v=B], scaled by 2*den:  sum_C (P1_C + P2_C) = 2 den X"""
    N, n = S['N'], S['n']; C = set(S['Cs'][ci]); w = int(S['al'][len(C)] * S['den'])
    mem = np.array([C <= set(A) for A in S['fam']])
    eAs = np.array(S['eAs'])                       # [A, i, r]
    Em = eAs * mem[:, None, None]
    # P1(A,i,B,j) = w sum_r e_A[i,r] e_B[j,r]  -> T[x=A,y=i,u=j,v=B]
    P1 = w * np.einsum('air,bjr->aijb', Em, Em)
    # P2(A,i,B,j) = w sum_r e_B[i,r] e_A[j,r] -> T[A,i,j,B]
    P2 = w * np.einsum('bir,ajr->aijb', Em, Em)
    return P1, P2
if __name__ == "__main__":
    r, k = int(sys.argv[1]), int(sys.argv[2]); order = sys.argv[3] if len(sys.argv) > 3 else 'lex'
    S = hr.setup(r, k, order=order); N, n = S['N'], S['n']
    X = np.zeros((N, n, n, N), dtype=np.int64)
    for x in range(N):
        for y in range(n): X[x, y, y, x] = 2 * S['den']
    P = [piece_tensors(S, ci) for ci in range(len(S['Cs']))]
    tot = sum(a + b for a, b in P); print("sum of pieces == 2den X:", (tot == X).all())
    T1 = sum(a for a, _ in P); T2 = sum(b for _, b in P)
    print(f"N={N} t={S['t']} k={k}: rr(T1)={rr(T1)} rc(T2)={rc(T2)} | rc(T1)={rc(T1)} rr(T2)={rr(T2)}")
    prof = []
    print(" C        M_C  rr(P1) rc(P1) | rr(P2) rc(P2) | rr(P1+P2) rc(P1+P2)")
    for ci, (a, b) in enumerate(P):
        C = S['Cs'][ci]; M = sum(set(C) <= set(A) for A in S['fam'])
        row = (C, M, rr(a), rc(a), rr(b), rc(b), rr(a + b), rc(a + b)); prof.append(row)
        print(f" {str(C):8} {M:4} {row[2]:6} {row[3]:6} | {row[4]:6} {row[5]:6} | {row[6]:8} {row[7]:8}", flush=True)
    # threshold: move C-side pieces P2 into R in order of increasing M_C (cheapest first)
    idx = sorted(range(len(P)), key=lambda c: prof[c][1])
    print("threshold (move P2_C to R, increasing M_C):  m  p  q")
    curve = []
    for m in range(len(P) + 1):
        moved = set(idx[:m])
        R = T1 + sum((P[c][1] for c in moved), np.zeros_like(X)); C2 = X - R
        curve.append((m, rr(R), rc(C2))); print("  ", curve[-1], flush=True)
    json.dump(dict(N=N, t=S['t'], k=k, order=order, prof=[list(map(lambda z: list(z) if isinstance(z, tuple) else z, p)) for p in prof], curve=curve),
              open(f"pieces_N{N}_k{k}_{order}.json", "w"))
