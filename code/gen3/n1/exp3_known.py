import numpy as np, math, itertools
from rho import *
print("eight-square identity (12): mismatches over 4096 settings =", verify_eight_square())
# eight-square pieces on one 3-bit block
N = 8
R = np.zeros((N,)*4); C = np.zeros((N,)*4)
for x, y, u, v in itertools.product(range(8), repeat=4):
    bx, by, bu, bv = B3[x], B3[y], B3[u], B3[v]
    if all((a ^ b ^ c ^ d) == 0 for a, b, c, d in zip(bx, by, bu, bv)):
        R[x, y, u, v] = sigma3(bx, by)*sigma3(bv, bu)
        C[x, y, u, v] = sigma3(bx, bu)*sigma3(bv, by)
assert np.array_equal(R + C, 2*crossing(8))
print("8-square: ranks R row/col", rank(rowmat(R)), rank(colmat(R)), " C row/col", rank(rowmat(C)), rank(colmat(C)))
# paper's simple count (Lemma 2): per block 3+6 bits + 1+1 naming bits = 11 of 12
print("paper count rho = 11/12 =", 11/12, "; Gibbs/entropy version log2(sum sqrt(r s))/(2 lam) =",
      math.log2(2*math.sqrt(8*64))/(2*3))
r, rc, cc, d, tau = rho_of(R, 2, 3); print(f"8-square, exact-rank (eigen) count: rho={r:.4f} (paper App.H: 10.81/12={10.81/12:.4f}) row{rc} col{cc} diag={d} tau*={tau}")
# four-square: restriction to blocks (p1,p2,0)
idx4 = [i for i in range(8) if B3[i][2] == 0]
R4 = R[np.ix_(idx4, idx4, idx4, idx4)]; C4 = C[np.ix_(idx4, idx4, idx4, idx4)]
assert np.array_equal(R4 + C4, 2*crossing(4))
r, rc, cc, d, tau = rho_of(R4, 2, 2); print(f"4-square: rho={r:.4f} (paper 7.58/8={7.58/8:.4f}) row{rc} col{cc}")
# single-bit identity of App. H: [x=v][y=u] = [x=u][y=v] - [x!=y](-1)^x [u!=v](-1)^u
R1 = np.zeros((2,)*4)
for x, y, u, v in itertools.product(range(2), repeat=4):
    R1[x, y, u, v] = (x == u)*(y == v)
C1 = crossing(2) - R1
r, rc, cc, d, tau = rho_of(R1, 1, 1); print(f"1-bit: rho={r:.4f} (paper log2(12)/4={math.log2(12)/4:.4f}) row{rc} col{cc}")
# validate eigen-count model against exact grouped ranks at finite t
for name, (A, Bm, lam) in {"1-bit": (R1, C1, 1), "4-sq": (R4, C4, 2)}.items():
    for t in ([2, 3, 4, 5] if lam == 1 else [2, 3]):
        gr = grouped_ranks(A, Bm, t)
        best = max(math.log2(a) + math.log2(b) for a, b in gr)
        print(f"  {name} t={t}: exact max_i log2(rank_row*rank_col)/(4 lam t) = {best/(4*lam*t):.4f}  ranks={gr}")
