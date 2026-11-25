# Which class groups suffice to realize Smith modes? rows from class set G (all columns): count new directions at s_j beyond global V(s_{j-1}).
import sys, itertools, math
from fractions import Fraction as Fr
from lib12 import *
from levelm3 import adapted_basis, K, MOD
from obstr import gf2rank
def run(a,t,zeta,groups):
    E,c,T,B=cauchy_block(a,t,zeta,0); n=len(E); i=v2(t)
    Ez=[[toZ(z,MOD) for z in r] for r in E]
    s0=a-i+1+v2q(zeta); s=[s0+2*j+2*v2(math.factorial(j)) for j in range(1<<(a-i-1))]
    G=adapted_basis(Ez)
    for grp in groups:
        rows=[k for k,x in enumerate(B) if x%(1<<i) in grp]
        A=adapted_basis([Ez[r] for r in rows])
        info=[]
        for j in range(len(s)):
            Vr=[red for v,red in A if v<=s[j]]; Vg=[red for v,red in G if v<=(s[j-1] if j else -1)]
            info.append(gf2rank(Vg+Vr)-gf2rank(Vg))
        print(a,t,"i",i,"group",grp,"rowSmith",[v for v,_ in A],"new dirs per mode j (want |group| each):",info)
run(5,4,1,[(0,),(0,2),(0,1),(0,1,2,3)])
run(6,8,1,[(0,4),(0,2,4,6),(0,1,2,3,4,5,6,7)])
run(6,4,1,[(0,2),(0,1,2,3)])
