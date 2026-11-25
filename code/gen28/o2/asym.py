# Asymptotic (k -> inf) cost of type-class pieces of the k-th tensor power of ONE piece (SWAP = R + C), from the spectra of
#   M_row = R_row S^-1 (eigen-0: n0, eigen-1: n1, other: ng)   and   M_col = C_col S^-1 (n0', n1', ng').
# Exact merged rank (proved for semisimple M with all 'other' eigenvalues outside (0,1) -- true for every piece tested):
#   rank_row(A_i) = #{words: a eigen-0 positions, b eigen-1, c other, a <= i <= a+c}   (i = number of C factors)
#   rank_col(A_i) = #{words: a' <= k-i <= a'+c'}  (eigen-0 of M_col forces R, eigen-1 forces C)
# h_merged(mu) = max_x [mu phi_row(x) + (1-mu) phi_col(x)] / ln N^2 ;  compared with
#   h_LB  = log_{N^2}(p^mu (N^2-q)^(1-mu) + (N^2-p)^mu q^(1-mu))   (lower bound for ANY type-homogeneous grouping)
#   h_O1  = log_{N^2}(p^mu N^(2(1-mu)) + N^(2 mu) q^(1-mu))          (G28-O1 entropy formula)
#   h_max = log_{N^2} max(p^mu N^(2(1-mu)), N^(2 mu) q^(1-mu))        (multiplicity-free, finite_hmax)
#   circle sqrt(mu^2 + (1-mu)^2)
import numpy as np, sys, json
from pieces import *
G = 400
al = np.linspace(0, 1, G+1)
A, C_ = np.meshgrid(al, al, indexing='ij'); B = 1 - A - C_; ok = B >= -1e-12; B = np.clip(B, 0, 1)
def xlnx(t): return np.where(t > 0, t*np.log(np.where(t > 0, t, 1)), 0.0)
def lnc(n): return np.log(n) if n > 0 else -np.inf
def phi_table(n0, n1, ng):
    H = -(xlnx(A) + xlnx(B) + xlnx(C_))
    with np.errstate(invalid='ignore'):
        V = H + np.where(A > 0, A*lnc(n0), 0) + np.where(B > 0, B*lnc(n1), 0) + np.where(C_ > 0, C_*lnc(ng), 0)
    V = np.where(ok, V, -np.inf); V = np.nan_to_num(V, nan=-np.inf)
    xs = al; out = np.empty(len(xs))
    for j, x in enumerate(xs):
        m = (A <= x + 1e-12) & (A + C_ >= x - 1e-12)
        out[j] = V[m].max() if m.any() else -np.inf
    return out
def spectrum_counts(M):
    ev = np.linalg.eigvals(M.astype(float)); n = len(ev)
    n0 = int((abs(ev) < 1e-7).sum()); n1 = int((abs(ev - 1) < 1e-7).sum())
    oth = ev[(abs(ev) >= 1e-7) & (abs(ev - 1) >= 1e-7)]
    inside = int(((oth.real > 1e-9) & (oth.real < 1 - 1e-9)).sum() + (abs(oth.imag) > 1e-7).sum())
    # semisimplicity of 0 and 1: rank(M) == rank(M^2), rank(M-1) == rank((M-1)^2)
    r = lambda X: np.linalg.matrix_rank(X.astype(float), tol=1e-7)
    I = np.eye(n); ss = (r(M) == r(M @ M)) and (r(M - I) == r((M - I) @ (M - I)))
    return n0, n1, n - n0 - n1, inside, ss, sorted(set(np.round(oth.real, 4)))
def compare(name, N, R, C, s, verbose=True):
    S = s*swap(N); N2 = N*N
    Mr = rowflat(R) @ np.linalg.inv(rowflat(S).astype(float)); Mc = colflat(C) @ np.linalg.inv(colflat(S).astype(float))
    cr = spectrum_counts(Mr); cc = spectrum_counts(Mc)
    p = N2 - cr[0]; q = N2 - cc[0]
    pr = phi_table(*cr[:3]); pc = phi_table(*cc[:3])[::-1]     # col side: constraint on 1-x
    mus = np.linspace(0, 1, 401); L = np.log(N2)
    hm = np.array([np.max(m*pr + (1-m)*pc) for m in mus])/L
    hLB = np.log(p**mus*(N2-q)**(1-mus) + (N2-p)**mus*q**(1-mus))/L
    hO1 = np.log(p**mus*N2**(1-mus) + N2**mus*q**(1-mus))/L
    hmax = np.log(np.maximum(p**mus*N2**(1-mus), N2**mus*q**(1-mus)))/L
    hB = np.hypot(mus, 1-mus)
    inn = slice(1, len(mus)-1); W = mus[inn]*(1-mus[inn])
    gi = lambda h: ((h-hB)[inn]/W).min()
    j = len(mus)//2
    if verbose:
        print(f"{name:28s} N={N:2d} (p,q)=({p},{q}) row(n0,n1,ng)={cr[:3]} col={cc[:3]} eig-in-(0,1)={cr[3]},{cc[3]} semisimple={cr[4]},{cc[4]}")
        print(f"     h(1/2): merged {hm[j]:.4f}  LB {hLB[j]:.4f}  O1-entropy {hO1[j]:.4f}  mult-free {hmax[j]:.4f}  circle 0.7071 |"
              f" min over interior mu of (h-circle)/(mu(1-mu)): merged {gi(hm):+.4f} LB {gi(hLB):+.4f} O1 {gi(hO1):+.4f} mult-free {gi(hmax):+.4f}")
    return mus, hm, hLB, hO1, hmax
if __name__ == '__main__':
    for which in sys.argv[1:]:
        if which == 'oct': N, (R, C, s) = octonion()
        elif which.startswith('F1_'): N, (R, C, s) = F1(int(which[3:]))
        else: N, (R, C, s) = load_json(which)
        compare(which.split('/')[-1], N, R, C, s)
