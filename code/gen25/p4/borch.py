# exact check of the Borcherds split on M_E0:  o(a_(-1)b) = sum_j a_(-1-j) b_(-1+j) + sum_j b_(-2-j) a_(j),  a=e^g1, b=e^g2
import sys, itertools
from modes import Setup, Fr
from objs import Copy
def mm(A, B):
    if not A or not B: return None
    return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def add(A, B):
    if B is None: return A
    return [[x+y for x, y in zip(r, s)] for r, s in zip(A, B)]
n, E0 = map(int, sys.argv[1:3])
St = Setup(n, E0); C = St.C
d = len(St.Mb[E0]); print("n", n, "E0", E0, "dims M_0..M_E0", [len(St.Mb[E]) for E in range(E0+1)], "S", St.S)
A = St.A; ia = {p: i for i, (p, m) in enumerate(A)}
Bs = {g: C.B(g, 0, 0) for g in St.S}
ok = True
for g1, g2 in itertools.product(St.S, repeat=2):
    i1 = ia[tuple(Fr(x) for x in g1)]; i2 = ia[tuple(Fr(x) for x in g2)]
    lhs = [[Fr(0)]*d for _ in range(d)]
    for g in St.S:
        c = Bs[g][i1][i2]
        if c:
            O = St.mode(g, E0, E0)
            lhs = [[x + c*y for x, y in zip(r, s)] for r, s in zip(lhs, O)]
    Pa = [[Fr(0)]*d for _ in range(d)]; Pb = [[Fr(0)]*d for _ in range(d)]
    for j in range(E0+1):
        Pa = add(Pa, mm(St.mode(g1, E0-j, E0), St.mode(g2, E0, E0-j)))
    for j in range(E0):
        Pb = add(Pb, mm(St.mode(g2, E0-1-j, E0), St.mode(g1, E0, E0-1-j)))
    rhs = add(Pa, Pb)
    if rhs != lhs: ok = False; print("MISMATCH", g1, g2)
print("Borcherds split exact on all pairs:", ok)
