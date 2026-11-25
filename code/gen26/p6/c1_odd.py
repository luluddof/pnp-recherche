# Direct check (no use of Prop. 3): full flattening F of one packet, a=3, i=0, zeta=1; R_p(F) at odd p.
import sys
src = open(sys.argv[1]).read().split("def cnt")[0]
exec(src)
a, i = 3, 0
L = 1 << a; t = 1 << i; T = t + (1 << (a + 1)); W = Wmat(a, T)
rg = [x for x in range(L) if (x >> i) & 1 == 0]; cg = [y for y in range(L) if (y >> i) & 1 == 1]
rh = [x for x in range(L) if (x >> i) & 1 == 1]; ch = [y for y in range(L) if (y >> i) & 1 == 0]
F = smith_vals(flattening(W, L, rg, cg, rh, ch), 30)
print("first invariants of F:", F[:20])
print("R_9(F) =", sum(1 for v in F if v < 9), "; Smith-monomial count of weight < 9 = 8 (1, Z0,Z0',Z1,Z1', Z0^2, Z0Z0', Z0'^2)")
