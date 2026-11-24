# rigorous upper bound on rank F^{(n)}: sum over charge blocks q in Z^n of min(rows_q, cols_q),
# rows_q = [x^E0 y^Ha] prod_i R_{q_i}(x,y), cols_q = [x^E0 y^Hb] prod_i C_{q_i}(x,y)   (per-copy weights NOT shared)
import itertools, math
from pred import part, dV, dM, CH, comps
def rowpoly(q,Emax,Hmax):  # dict (e,ha) -> rows for one copy, charge block q
    return {(e,h): sum(dV(h,s)*dM(e,q+s) for s in CH) for e in range(Emax+1) for h in range(Hmax+1)}
def colpoly(q,Emax,Hmax):
    return {(e,h): sum(dV(h,t)*dM(e,q-t) for t in CH) for e in range(Emax+1) for h in range(Hmax+1)}
def coef(polys,E0,H):
    # coefficient of x^E0 y^H in product of 2-var polys
    cur={(0,0):1}
    for P in polys:
        new={}
        for (e1,h1),c1 in cur.items():
            for (e2,h2),c2 in P.items():
                if c2 and e1+e2<=E0 and h1+h2<=H: new[(e1+e2,h1+h2)]=new.get((e1+e2,h1+h2),0)+c1*c2
        cur=new
    return cur.get((E0,H),0)
def bound(n,E0,ha,hb):
    tot=0
    for q in itertools.product(CH,repeat=n):
        X=coef([rowpoly(qq,E0,ha) for qq in q],E0,ha); Y=coef([colpoly(qq,E0,hb) for qq in q],E0,hb)
        tot+=min(X,Y)
    return tot
if __name__=="__main__":
    exact={(0,0,0):4,(1,1,1):24,(2,0,0):52,(0,2,2):48,(2,1,1):104,(2,2,2):408,(1,2,2):104,(2,1,3):198,(3,1,1):208,(3,2,2):820,(2,3,3):818,(3,3,3):1640,(4,2,2):1866,(4,1,1):480,(4,3,3):3730,(4,4,4):8470}
    for k,v in exact.items(): print("2copy",k,"exact",v,"bound",bound(2,*k))
