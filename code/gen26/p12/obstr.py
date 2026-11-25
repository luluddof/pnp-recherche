# Obstruction check: at level i (D=2^{a-i}), compare Smith of class-rho rows (all columns) with class profile s_j, and test
# whether V_rho(s_j) (reductions of class-rho-row elements with val<=s_j) is contained in global V(s_{j-1}).
import sys, itertools, math
from fractions import Fraction as Fr
from lib12 import *
from levelm3 import adapted_basis, K, MOD
def gf2rank(vs):
    bas={}
    for x in vs:
        while x:
            h=x.bit_length()-1
            if h in bas: x^=bas[h]
            else: bas[h]=x; break
    return len(bas)
def run(a,t,zeta,b):
    E,c,T,B=cauchy_block(a,t,zeta,b); n=len(E); i=v2(t)
    Ez=[[toZ(z,MOD) for z in r] for r in E]
    s0=a-i+1+v2q(zeta); s=[s0+2*j+2*v2(math.factorial(j)) for j in range(1<<(a-i-1))]
    G=adapted_basis(Ez)
    out=[]
    for rho in sorted(set(x%(1<<i) for x in B)) if i>0 else [0]:
        rows=[k for k,x in enumerate(B) if (x%(1<<i)==rho if i>0 else True)]
        A=adapted_basis([Ez[r] for r in rows])
        vals=[v for v,_ in A]
        # dim of class-local reductions at <= s_j inside global V(s_{j-1})
        info=[]
        for j in range(1,len(s)):
            Vr=[red for v,red in A if v<=s[j]]; Vg=[red for v,red in G if v<=s[j-1]]
            new=gf2rank(Vg+Vr)-gf2rank(Vg)
            info.append(new)
        out.append((rho,vals,info))
    print(a,t,zeta,"i",i,"s",s,"global",[v for v,_ in G][:12],"...")
    for rho,vals,info in out[:2]: print("   class",rho,"rowSmith",vals,"new dirs beyond V(s_{j-1}) at s_j:",info)
for (a,t,z) in [(4,2,1),(4,6,1),(4,2,3),(4,2,2),(5,4,1),(5,2,1),(5,1,1),(4,1,1),(5,8,1),(6,8,1),(6,4,1)]:
    run(a,t,Fr(z),0)
