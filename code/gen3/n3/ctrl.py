import sys; sys.path.insert(0,'.')
from nsz import *
# Tseitin on K4 (6 edges, vertex degree 3), odd total charge -> UNSAT
E=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
cl=tseitin(E,{0:1}); n=6
print('tseitin K4', degs(n,[clause_poly(c) for c in cl],6))
# Tseitin on 3x3 grid-ish: Petersen-like small 3-regular graph K_{3,3}
E=[(a,b) for a in range(3) for b in range(3,6)]
cl=tseitin(E,{0:1}); n=9
print('tseitin K33', degs(n,[clause_poly(c) for c in cl],9))
