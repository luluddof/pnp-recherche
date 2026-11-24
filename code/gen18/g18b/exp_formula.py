#!/usr/bin/env python3
"""G18b exp_formula: end-to-end check at every node t of the mod-3 evaluator on a tiny 3-CNF (<=3 occurrences),
coefficients in F_9, per-node measurement matrices C_{t,g} from the piecewise construction (r=2, Q=4, m=3).
Checks: (a) each replacement form is a function of (x_g, y_{-g} mod 3); (b) sum over reachable y of prod_g B_g(y)
equals brute-force H_w(t) (weighted, weights in F_9); (c) H_w(0) = weighted count; (d) #tuples <= 3^{sum rank_3}."""
import itertools, math, random, sys
import numpy as np
from exp_aux import *
def main(seed):
    random.seed(seed)
    Q,m=4,3
    nodes=[f9(1,0),f9(2,0),f9(0,1),f9(0,2),f9(1,1),f9(1,2),f9(2,1),f9(2,2)]
    obj=build(Q,m,nodes); D=obj['D']; signs=obj['signs']; Psi=obj['Psi']
    n=10
    # random 3-CNF, each variable in at most 3 clauses
    occ=[v for v in range(n) for _ in range(3)]; random.shuffle(occ)
    clauses=[]
    while len(occ)>=3 and len(clauses)<9:
        vs=[]
        for v in list(occ):
            if v not in vs: vs.append(v); occ.remove(v)
            if len(vs)==3: break
        if len(vs)<3: break
        clauses.append([(v,random.random()<0.5) for v in vs])
    a=[[random.randrange(2) for _ in range(Q)] for _ in range(n)]
    d=[[random.randrange(2) for _ in range(Q)] for _ in clauses]
    lab=lambda e,v: sum(4**s*(2*a[v][s]+d[e][s]) for s in range(Q))
    place={}; used={}
    for e,cl in enumerate(clauses):
        for v,_ in cl:
            u=lab(e,v); c=used.get(u,0); used[u]=c+1; place[(e,v)]=(c,u)
    h=max(used.values())
    w=[nodes[random.randrange(8)] for _ in range(n)]            # random nonzero weights in F_9
    X=list(itertools.product((0,1),repeat=n))
    def wt(x):
        r=ONE
        for i in range(n):
            if x[i]: r=mul(r,w[i])
        return r
    def phi(vals,cl):   # 1 - prod (1 - lit)
        pr=ONE
        for val,(v,neg) in zip(vals,cl):
            lit=add(ONE,mul(f9(2,0),val)) if neg else val
            pr=mul(pr,add(ONE,mul(f9(2,0),lit)))
        return add(ONE,mul(f9(2,0),pr))
    def Hbrute(Pt):
        # rows over variables
        rows={}
        for (e,v),(c,u) in place.items():
            row=[ZERO]*n
            for (e2,v2),(c2,u2) in place.items():
                if c2==c: row[v2]=add(row[v2],(int(Pt[0][u,u2]),int(Pt[1][u,u2])))
            rows[(e,v)]=row
        tot=ZERO; vals_all={}
        for x in X:
            pr=wt(x); vv={}
            for e,cl in enumerate(clauses):
                vals=[]
                for v,_ in cl:
                    s=ZERO
                    for i in range(n):
                        if x[i]: s=add(s,rows[(e,v)][i])
                    vals.append(s)
                vv[e]=tuple(vals); pr=mul(pr,phi(vals,cl))
            vals_all[x]=vv; tot=add(tot,pr)
        return tot, vals_all
    I=np.eye(D,dtype=np.int64)
    H0,_=Hbrute([I,np.zeros((D,D),dtype=np.int64)])
    cnt=ZERO; nsat=0
    for x in X:
        if all(any((x[v]==0) if neg else (x[v]==1) for v,neg in cl) for cl in clauses): cnt=add(cnt,wt(x)); nsat+=1
    assert H0==cnt
    print(f"seed {seed}: n={n}, {len(clauses)} clauses, h={h}, #sat={nsat}; H_w(0) = weighted count in F_9: OK")
    for t in range(Q):
        tau=nodes[t]; Pt=obj['Pmat'](tau)
        Ht,vals_all=Hbrute(Pt)
        badset={z for z in set(signs) if add(peval(obj['Pz'][z],tau),
                 peval(obj['Pz'][tuple(-z[j] if j==t else z[j] for j in range(Q))],tau))!=ZERO}
        brow=[k for k,z in enumerate(signs) if z in badset]
        grp=[a[v][t] for v in range(n)]
        Cg={0:[],1:[]}
        for c in range(h):
            for k in brow:
                row=np.zeros(n,dtype=np.int64)
                for (e2,v2),(c2,u2) in place.items():
                    if c2==c: row[v2]+=Psi[u2,k]
                for g in (0,1):
                    rg=np.where(np.array(grp)==g,row,0)
                    if np.abs(rg).any(): Cg[g].append(rg)
        Cm={g:(np.array(Cg[g]) if Cg[g] else np.zeros((0,n),dtype=np.int64)) for g in (0,1)}
        rk={g:rank_mod3(Cm[g])[0] if len(Cm[g]) else 0 for g in (0,1)}
        yof=lambda x,g: tuple((Cm[g]@np.array(x))%3) if len(Cm[g]) else ()
        yofZ=lambda x,g: tuple((Cm[g]@np.array(x))) if len(Cm[g]) else ()
        xg=lambda x,g: tuple(x[i] for i in range(n) if grp[i]==g)
        # (a) forms are functions of (x_g, y_other)
        table={}; ok=True
        for x in X:
            ys=(yof(x,0),yof(x,1))
            for e in range(len(clauses)):
                g=d[e][t]; key=(e,xg(x,g),ys[1-g])
                if table.setdefault(key,vals_all[x][e])!=vals_all[x][e]: ok=False
        # (b) evaluator
        Y={g:sorted({yof(x,g) for x in X}) for g in (0,1)}
        YZ={g:len({yofZ(x,g) for x in X}) for g in (0,1)}
        U={g:sorted({xg(x,g) for x in X}) for g in (0,1)}
        gi={g:[i for i in range(n) if grp[i]==g] for g in (0,1)}
        yg_of_u={g:{u:tuple((Cm[g][:,gi[g]]@np.array(u,dtype=np.int64))%3) if len(Cm[g]) and len(u) else
                    (tuple([0]*len(Cm[g]))) for u in U[g]} for g in (0,1)}
        Z=ZERO
        for y0 in Y[0]:
            for y1 in Y[1]:
                ys=(y0,y1); pr=ONE
                for g in (0,1):
                    Bg=ZERO
                    for u in U[g]:
                        if yg_of_u[g][u]!=ys[g]: continue
                        aa=ONE
                        for j,i in enumerate(gi[g]):
                            if u[j]: aa=mul(aa,w[i])
                        for e,cl in enumerate(clauses):
                            if d[e][t]==g: aa=mul(aa,phi(table[(e,u,ys[1-g])],cl))
                        Bg=add(Bg,aa)
                    pr=mul(pr,Bg)
                Z=add(Z,pr)
        print(f"  t={t}: rows per group {[len(Cm[0]),len(Cm[1])]} rank_3 {[rk[0],rk[1]]} | tuples mod 3: {len(Y[0])*len(Y[1])}"
              f" <= 3^(sum rank) = {3**(rk[0]+rk[1])}; reachable over Z: {YZ[0]*YZ[1]}; group sizes {[len(gi[0]),len(gi[1])]}"
              f" | forms depend only on (x_g, y_-g mod 3): {'OK' if ok else 'FAIL'} | evaluator == brute H_w(t): {'OK' if Z==Ht else 'FAIL'}")
        assert ok and Z==Ht and len(Y[0])*len(Y[1])<=3**(rk[0]+rk[1])
if __name__=="__main__":
    for s in (1,2,3): main(s)
