# Johnson-regraded Hrubes (JH) family: X = T1 + T2, T1 = F1(A,B) e_A e_B^t, T2 = F2(A,B) e_B e_A^t (swap), F1 + E o F2 = I.
# p = N rank F1, q = N rank F2, ||Fl1 R|| = ||F1||, ||Fl2 C|| = ||F2|| (checked on full tensors for N <= 32).
import sys, itertools, json, math, numpy as np
from fractions import Fraction
import johnson
from cliff import rho
from spread import *
def FF(t, k, N, j1, j2):
    fam = list(itertools.combinations(range(t), k))[:N]
    phi1, phi2 = johnson.solve(k, j1, j2)
    M = np.array([[len(set(a) & set(b)) for b in fam] for a in fam])
    f1 = np.array([float(x) for x in phi1]); f2 = np.array([float(x) for x in phi2])
    return fam, f1[M], f2[M], M
def full_check(r, k, j1, j2):
    import hr
    S = hr.setup(r, k); N = S['N']; eAs = np.array(S['eAs'], dtype=float)
    fam, F1, F2, M = FF(S['t'], k, N, j1, j2)
    T2 = np.einsum('ab,bir,ajr->aijb', F2, eAs, eAs)
    R = swapT(N) - T2
    d = analyze(R, T2)
    return d, np.linalg.norm(F1, 2), np.linalg.norm(F2, 2)
if __name__ == "__main__":
    if sys.argv[1] == 'check':
        for r, k in [(3, 2), (4, 2), (5, 2)]:
            for j1 in range(k):
                d, n1, n2 = full_check(r, k, j1, k - 1 - j1)
                print(f"r={r} k={k} j=({j1},{k-1-j1})", fmt(d), f" ||F1||={n1:.4g} ||F2||={n2:.4g}", flush=True)
        sys.exit()
    out = []
    for r in range(3, 13):
        N = 2 ** r; t = rho(N) - 1
        ks = [k for k in range(2, 7) if math.comb(t, k) >= N][:2]
        for k in ks:
            for j1 in range(k):
                j2 = k - 1 - j1
                fam, F1, F2, M = FF(t, k, N, j1, j2)
                assert np.allclose(F1 + (-1.0) ** (k - M) * F2, np.eye(N))
                e1 = np.linalg.eigvalsh(F1); e2 = np.linalg.eigvalsh(F2)
                r1 = int((np.abs(e1) > 1e-8 * np.abs(e1).max()).sum()); r2 = int((np.abs(e2) > 1e-8 * np.abs(e2).max()).sum())
                p, q = N * r1, N * r2; o1, o2 = np.abs(e1).max(), np.abs(e2).max()
                L = max(o1 / math.sqrt(q), o2 / math.sqrt(p)); rat = N * N / (p * o1 + q * o2)
                rec = dict(N=N, t=t, k=k, j1=j1, j2=j2, p=p, q=q, th=math.log(p, N), T=math.log(q, N), opR=o1, opC=o2, L=L, ineq=rat,
                           a1=math.log(o1 / math.sqrt(q), N), a2=math.log(o2 / math.sqrt(p), N))
                out.append(rec)
                print(f"N={N:5d} t={t} k={k} j=({j1},{j2}) p={p:9d} q={q:9d} (th,T)=({rec['th']:.3f},{rec['T']:.3f}) ||R||={o1:9.4g} ||C||={o2:9.4g} L={L:.4g} "
                      f"log_N(|R|/sqrt q)={rec['a1']:.3f} log_N(|C|/sqrt p)={rec['a2']:.3f} N2/(p|R|+q|C|)={rat:.4f}", flush=True)
    json.dump(out, open('jh.json', 'w'), indent=0)
