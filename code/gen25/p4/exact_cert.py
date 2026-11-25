# EXACT counterexample to "the P_b families are always triangularisable":
#  n copies, N = n, a_x = e^{alpha_x}, b_v = e^{sum_{i != v} alpha_i}  (so a_x(-1) b_v = delta_xv e^{alpha_[n]} up to sign, B_g exact),
#  labels: Y = random integer columns in M_E0, U = exact left inverse of (+-O_{[n]}) Y  => identity holds exactly.
#  Certificate: dim_Q span{E_m} > N(N+1)/2 (and same for C_m) => not triangularisable by equivalence over any extension of Q
#  (an upper-triangular-after-equivalence family spans a space of dim <= N(N+1)/2).
import sys, random, itertools
import flint
from modes import Setup, Fr
n, E0, seed = map(int, sys.argv[1:4])
random.seed(seed)
St = Setup(n, E0); C = St.C; N = n
def fq(x): return flint.fmpq(x.numerator, x.denominator)
def Q(M): return flint.fmpq_mat([[fq(x) if isinstance(x, Fr) else flint.fmpq(x) for x in r] for r in M])
A = St.A; ia = {tuple(int(c) for c in p): i for i, (p, m) in enumerate(A)}
full = tuple([1]*n)
ag = [tuple(1 if i == x else 0 for i in range(n)) for x in range(N)]
bg = [tuple(0 if i == v else 1 for i in range(n)) for v in range(N)]
Bs = {g: C.B(g, 0, 0) for g in St.S}
beta = {(g, x, v): Bs[g][ia[ag[x]]][ia[bg[v]]] for g in St.S for x in range(N) for v in range(N)}
nz = {(g, x, v): b for (g, x, v), b in beta.items() if b}
print("nonzero pair coefficients:", nz)
O = Q(St.mode(full, E0, E0)); dM = O.nrows()
print("n=%d E0=%d dim M_E0=%d rank O_[n] on M_E0 = %d" % (n, E0, dM, O.rank()))
if O.rank() < N: print("rank too small for N=n"); sys.exit()
sgn = fq(list(nz.values())[0])
Y = flint.fmpq_mat([[random.randint(-3, 3) for _ in range(N)] for _ in range(dM)])
OY = O*Y*sgn
G = OY.transpose()*OY
U = G.inv()*OY.transpose()      # N x dM, U O Y sgn = I
# identity check, exact
ok = True
Os = {g: Q(St.mode(g, E0, E0)) for g in St.S}
for x, v in itertools.product(range(N), repeat=2):
    tot = flint.fmpq_mat(N, N)
    for g in St.S:
        b = beta[(g, x, v)]
        if b: tot += U*Os[g]*Y*fq(b)
    want = flint.fmpq_mat(N, N, [1 if (x == v and i == j) else 0 for i in range(N) for j in range(N)])
    ok &= (tot == want)
print("identity sum_g beta_g(x,v) U O_g Y = delta_xv I exact:", ok)
# families: E_m(x,y) = <m*, (a_x)_(j) Y_y>,  C_m(v,u) = <U_u, (b_v)_(-2-j) m>,  j = E0-1-wt(m)
Erows = []; Crows = []
for E in range(E0):
    Dn = {x: Q(St.mode(ag[x], E0, E))*Y for x in range(N)}      # dim M_E x N
    Up = {v: U*Q(St.mode(bg[v], E, E0)) for v in range(N)}      # N x dim M_E
    for m in range(len(St.Mb[E])):
        Erows.append([Dn[x][m, y] for x in range(N) for y in range(N)])
        Crows.append([Up[v][u, m] for v in range(N) for u in range(N)])
dE = flint.fmpq_mat(Erows).rank(); dC = flint.fmpq_mat(Crows).rank()
print("dim M_<E0 = %d ; dim span{E_m} = %d, dim span{C_m} = %d ; bound for triangularisable N(N+1)/2 = %d" % (len(Erows), dE, dC, N*(N+1)//2))
print("VERDICT: E-family %s, C-family %s" % ("NOT triangularisable" if dE > N*(N+1)//2 else "inconclusive by dim", "NOT triangularisable" if dC > N*(N+1)//2 else "inconclusive by dim"))
