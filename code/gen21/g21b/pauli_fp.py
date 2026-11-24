# Pauli-diagonal sub-family (G20c section 3) over GF(P), P odd, N = 2^n: (p,q) feasible in the family iff
# chi = g + f with |supp g| <= p and f = sum_{b in B} c_b omega_b, |B| <= q  (omega_b(a) = (-1)^{omega(a,b)} symplectic).
# Exhaustive over B and c in GF(P)^q:  least p for each q.  (Over Q: diag.py of G19b gives 6,6,6,4,4 for q=1..5, N=4.)
# Dictionary valid for P odd: rank Phi_g = |supp g|, osr Phi_f = |supp f^| (needs 1/N^2 and W_b (x) W_b^-1 independent).
import sys, itertools, numpy as np
def run(n, P, qmax):
    pts = list(itertools.product([0, 1], repeat=2*n))
    def symp(a, b): return sum(a[i]*b[n+i] + a[n+i]*b[i] for i in range(n)) % 2
    chi = np.array([(-1)**(sum(a[i]*a[n+i] for i in range(n)) % 2) % P for a in pts])
    Om = np.array([[(-1)**symp(a, b) % P for a in pts] for b in pts])  # rows b
    M = len(pts); best = {}
    for q in range(1, qmax+1):
        cs = np.array(list(itertools.product(range(1, P), repeat=q)))  # nonzero coefficients (q exact; monotone anyway)
        bp = M
        for B in itertools.combinations(range(M), q):
            F = (cs @ Om[list(B)]) % P                 # all f for this B
            mism = (F != chi).sum(axis=1).min()
            if mism < bp: bp = mism
        best[q] = int(bp); print(f"n={n} N={2**n} GF({P}) q={q}: least p = {bp}", flush=True)
    return best
if __name__ == "__main__":
    n, P, qmax = map(int, sys.argv[1:4]); run(n, P, qmax)
