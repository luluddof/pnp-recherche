# After sparsify.py: threshold small entries, use the diagonal (torus) part of the gauge group to bring all nonzero entries to
# modulus 1 (linear least squares in logarithms), round to small rationals, verify EXACTLY with verify.py.
import numpy as np, sys
from fractions import Fraction
from verify import verify
N, p, q = map(int, sys.argv[1].split(',')); thr = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
z = np.load(f'sparse_{N}_{p}_{q}.npy'); G = z[:q*N*N].reshape(q,N,N).copy(); Gp = z[q*N*N:].reshape(q,N,N).copy()
G[np.abs(G) < thr] = 0; Gp[np.abs(Gp) < thr] = 0
rows = []; rhs = []      # unknowns: mu_l (q), gam_x (N), eta_u (N)
for l in range(q):
    for a in range(N):
        for b in range(N):
            if G[l,a,b] != 0:
                r = np.zeros(q+2*N); r[l] = 1; r[q+a] += 1; r[q+N+b] -= 1; rows.append(r); rhs.append(-np.log(abs(G[l,a,b])))
            if Gp[l,a,b] != 0:
                r = np.zeros(q+2*N); r[l] = -1; r[q+N+a] += 1; r[q+b] -= 1; rows.append(r); rhs.append(-np.log(abs(Gp[l,a,b])))
sol = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)[0]
mu, gam, eta = np.exp(sol[:q]), np.exp(sol[q:q+N]), np.exp(sol[q+N:])
A = np.einsum('l,x,lxu,u->lxu', mu, gam, G, 1/eta); B = np.einsum('l,y,lyv,v->lyv', 1/mu, eta, Gp, 1/gam)
fr = lambda t: Fraction(float(t)).limit_denominator(12)
Af = [[[fr(A[l,a,b]) for b in range(N)] for a in range(N)] for l in range(q)]; Bf = [[[fr(B[l,a,b]) for b in range(N)] for a in range(N)] for l in range(q)]
err = max(max(abs(float(Af[l][a][b]) - A[l,a,b]), abs(float(Bf[l][a][b]) - B[l,a,b])) for l in range(q) for a in range(N) for b in range(N))
print("max rounding error", err, " distinct values:", sorted(set(str(e) for M in (Af, Bf) for l in M for r in l for e in r)))
pp, qq = verify(N, Af, Bf, f"N={N} target ({p},{q})")
def show(M): return [" ".join(f"{'+' if e > 0 else '-'}{'' if abs(e) == 1 else abs(e)}E{a+1}{b+1}" for a in range(N) for b in range(N) if (e := M[l][a][b]) != 0) for l in range(q)]
print("G  :", show(Af)); print("G' :", show(Bf))
import json; json.dump({"N": N, "p": pp, "q": qq, "G": [[[str(e) for e in r] for r in l] for l in Af], "Gp": [[[str(e) for e in r] for r in l] for l in Bf]}, open(f"exact_{N}_{pp}_{qq}.json", "w"))
