import sys; sys.path.insert(0,'/tmp/claude-0/-home-user/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/step3d')
from ns import *
c=0
for s in range(200):
    cl=rand3cnf(8,36,s)
    if unsat(8,cl): continue
    P=[clause_poly(x) for x in cl]
    for q in (2,3,4): assert min_deg(8,P,q,8) is None
    c+=1
    if c==10: break
print('sat instances never refuted:',c)
