import sys; sys.path.insert(0,'/tmp/claude-0/-home-user/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/step3d')
from ns import *
K4=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
K33=[(a,b) for a in range(3) for b in range(3,6)]
Q3=[(a,a^(1<<k)) for a in range(8) for k in range(3) if a<a^(1<<k)]
for name,E in [('K4',K4),('K3,3',K33),('Q3',Q3)]:
    cl=tseitin(E,{0:1}); n=len(E); assert unsat(n,cl)
    P=[clause_poly(c) for c in cl]
    print(name,'n=',n,{q:min_deg(n,P,q,n) for q in (2,3,4)},flush=True)
