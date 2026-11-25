# find numerical true-N solutions with G24f's numsearch, build the P_b families E_m, C_m, and test triangularisability
import sys, numpy as np, pickle, time
sys.path.insert(0, '/home/user/pnp-recherche/code/gen24/g24f')
import numsearch
from modes import Setup, npf
n, E0, N, trials, seed = map(int, sys.argv[1:6])
out = sys.argv[6]
St = Setup(n, E0)
Bn, On, dA, dB, dM = numsearch.prep(n, E0, 0, 0)
# numsearch.prep drops channels with zero B or O; we need A in the C.basis(0) coordinates (same order as St.A)
sols = []
for s in range(seed, seed + trials):
    t = time.time()
    best, z, unpack = numsearch.solve(Bn, On, dA, dB, dM, N, trials=1, seed=s, tol=1e-9)
    A, Bm, U, Y = unpack(z)
    print("seed", s, "res %.2e maxabs %.2f (%.0fs)" % (best, np.abs(z).max(), time.time()-t), flush=True)
    if best < 1e-7: sols.append(dict(A=A, Bm=Bm, U=U, Y=Y, res=best))
# families
gs = [tuple(int(x) for x in p) for (p, m) in St.A]
modes_down = {(g, E): npf(St.mode(g, E0, E)) for g in gs for E in range(E0)}   # e^g_(E0-1-E): M_E0 -> M_E
modes_up = {(g, E): npf(St.mode(g, E, E0)) for g in gs for E in range(E0)}     # e^g_(-1-E0+E): M_E -> M_E0
res = []
for so in sols:
    A, Bm, U, Y = so['A'], so['Bm'], so['U'], so['Y']
    Es = []; Cs = []; lab = []
    for E in range(E0):
        # E_m(x,y) = sum_g A[x,g] (e^g_(j) Y)[m,y] ; C_m(v,u) = sum_g Bm[v,g] (U e^g_(-2-j))[u,m]
        Tm = np.einsum('xg,gmy->mxy', A, np.array([modes_down[(g, E)] @ Y for g in gs]))
        Cm = np.einsum('vg,gum->mvu', Bm, np.array([U @ modes_up[(g, E)] for g in gs]))
        for k in range(Tm.shape[0]): Es.append(Tm[k]); Cs.append(Cm[k]); lab.append(E)
    Es = np.array(Es); Cs = np.array(Cs)
    Pb = np.einsum('mxy,mvu->xvuy', Es, Cs)
    # full o(a_x(-1) b_v) sandwiched
    beta = np.einsum('xi,gij,vj->gxv', A, Bn, Bm); Om = np.einsum('uk,gkl,ly->guy', U, On, Y)
    T = np.einsum('gxv,guy->xvuy', beta, Om)
    I = np.einsum('xv,uy->xvuy', np.eye(N), np.eye(N))
    Pa = T - Pb
    sv = np.linalg.svd(Pa.transpose(0, 2, 1, 3).reshape(N*N, N*N), compute_uv=False)  # (x,u)|(v,y)
    print("  |T-I| %.1e  rank_(xu|vy)(P_a) ~ %d (sv %s)" % (np.abs(T-I).max(), (sv > 1e-6).sum(), np.round(sv, 3)))
    res.append(dict(Es=Es, Cs=Cs, lab=lab, sol=so))
pickle.dump(res, open(out, 'wb'))
