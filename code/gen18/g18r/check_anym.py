# (#) is solvable over Z for EVERY m (Newton forward differences): so m need not be a power of 3. Scalar check over GF(27).
import math
from gf import GF
from design import Design
F=GF(3)
for (Q,m) in [(10,2),(10,4),(10,5),(11,6),(12,7)]:
    res=[]
    for shift in range(1,m+1):
        des=Design(F,Q,m,shift,list(range(1,Q+1)))
        maxdeg,nbp,acc=des.scalar_checks()
        res.append((maxdeg,round(float(des.bad_fraction(2)),4),round(float(des.bad_fraction(3)),4)))
    print(f"Q={Q} m={m}: all shifts pass (P_z(0)=1, deg<=m-1, pairing inside blocks). (deg, frac r=2, frac r=3) per shift: {res}; 1/m={1/m:.4f}")
# integer solution: a_i = sum_k (-1)^(i-k) C(i,k) (-1)^(lo+k) composed with Vandermonde shift; check unimodularity of [C(lo+j,i)]
import sympy
for m in (4,5,6,7,10):
    for lo in (0,1,5,17):
        M=sympy.Matrix(m,m,lambda j,i: math.comb(lo+j,i)); assert M.det()==1
print("det [C(lo+j,i)]_{j,i<m} = 1 for m in (4,5,6,7,10), lo in (0,1,5,17): (#) solvable over Z, hence over every field")
