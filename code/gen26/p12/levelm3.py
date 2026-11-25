# Level a-3 (D=8): test Lemma W''' criterion per class rho:
#  e_{rho,j} = reductions of class-rho modes (flag of V_{rho,1}); for each T subset {0,1,2,3}: is wedge_{j in T} e_{rho,j}
#  in V_{rho,|T|}(sum_{j in T} s_j) + span{e-monomials with smaller total mode index}?
import sys, itertools, math
from fractions import Fraction as Fr
from lib12 import *
K=160; MOD=1<<K
def vv(z): z%=MOD; return K if z==0 else v2(z)
def perZ(Mx):
    m=len(Mx); dp={0:1}
    for r in range(m):
        nd={}
        for mask,val in dp.items():
            for j in range(m):
                if not mask>>j&1:
                    nm=mask|1<<j; nd[nm]=(nd.get(nm,0)+val*Mx[r][j])%MOD
        dp=nd
    return dp[(1<<m)-1]
def adapted_basis(rows):
    """rows: list of lists mod MOD. returns list of (v, reduction as dict col->1 bitset int over column indices)"""
    R=[r[:] for r in rows]; out=[]
    alive=list(range(len(R)))
    while alive:
        best=(K+1,None,None)
        for a in alive:
            for c,z in enumerate(R[a]):
                if z:
                    vz=v2(z)
                    if vz<best[0]: best=(vz,a,c)
        v,a,c=best
        if a is None: break
        alive.remove(a)
        piv=R[a][c]; u=piv>>v; ui=pow(u,-1,MOD)
        R[a]=[(z*ui)%MOD for z in R[a]]
        for b in alive:
            if R[b][c]:
                f=R[b][c]>>v
                R[b]=[(zb-f*za)%MOD for zb,za in zip(R[b],R[a])]
        red=0
        for cc,z in enumerate(R[a]):
            if (z>>v)&1: red|=1<<cc
        out.append((v,red))
    return out
def run(a,t,zeta,b,rho):
    E,c,T,B=cauchy_block(a,t,zeta,b); n=len(E); i=v2(t); assert a-i==3
    Ez=[[toZ(z,MOD) for z in r] for r in E]
    s0=4+v2q(zeta); s=[s0,s0+2,s0+6,s0+8]
    rows=[k for k,x in enumerate(B) if x%(1<<i)==rho]
    rows.sort(key=lambda k:B[k])
    assert len(rows)==4
    # global e-basis: modes of ALL classes (needed for coordinates)
    classes=sorted(set(x%(1<<i) for x in B))
    emodes={}  # (class, j) -> linear form bitset over n cols
    for cl in classes:
        rr=[k for k,x in enumerate(B) if x%(1<<i)==cl]
        ab=adapted_basis([Ez[r] for r in rr])
        vs=[v for v,_ in ab]
        if vs!=s: print("class",cl,"mode valuations",vs,"!= s",s); return None
        for j,(v,red) in enumerate(ab): emodes[(cl,j)]=red
    keys=sorted(emodes); idxsum={k:k[1] for k in keys}
    # check basis
    mat=[emodes[k] for k in keys]
    # GF2 inverse: express X_y in e basis: solve
    basis={}
    piv_rows=[]
    import copy
    # Build matrix Emat (rows e_k over cols y). Need B with X_y = sum_k B[y][k] e_k.
    # Gaussian elimination on augmented [Emat | I]
    aug=[(mat[k], 1<<k) for k in range(n)]
    for col in range(n):
        p=next((r for r in range(col,n) if aug[r][0]>>col&1),None)
        if p is None: print("e not a basis"); return None
        aug[col],aug[p]=aug[p],aug[col]
        for r in range(n):
            if r!=col and aug[r][0]>>col&1:
                aug[r]=(aug[r][0]^aug[col][0], aug[r][1]^aug[col][1])
    # now row col has e-combination giving X_col: X_col = sum_{k in aug[col][1]} e_k
    Bm=[aug[y][1] for y in range(n)]  # bitset over keys indices
    res={}
    for d in range(1,5):
        Rsubs=list(itertools.combinations(rows,d)); Csubs=list(itertools.combinations(range(n),d))
        P=[[perZ([[Ez[r][u] for u in U] for r in R]) for U in Csubs] for R in Rsubs]
        ab=adapted_basis(P)
        for Tset in itertools.combinations(range(4),d):
            w=sum(s[j] for j in Tset)
            V=[red for v,red in ab if v<=w]
            # e-coordinates of each V vector restricted to monomials with idxsum >= sigma
            sig=sum(Tset)
            # monomials S (tuples of key indices) with idxsum>=sig, size d
            def ecoords(red):
                co={}
                for ci,U in enumerate(Csubs):
                    if red>>ci&1:
                        # X_U = wedge_{y in U} (sum_k B[y][k] e_k) -> sum_S det(B[U,S]) e_S
                        lists=[[k for k in range(n) if Bm[y]>>k&1] for y in U]
                        for combo in itertools.product(*lists):
                            if len(set(combo))<d: continue
                            S=tuple(sorted(combo))
                            co[S]=co.get(S,0)^1
                return {S:1 for S,vv_ in co.items() if vv_ and sum(keys[k][1] for k in S)>=sig}
            target=tuple(sorted(keys.index((rho,j)) for j in Tset))
            vecs=[ecoords(r) for r in V]
            # GF2: is unit(target) in span(vecs)?
            allS=sorted(set(S for vv_ in vecs for S in vv_)|{target}); ix={S:q for q,S in enumerate(allS)}
            bits=[sum(1<<ix[S] for S in vv_) for vv_ in vecs]
            tb=1<<ix[target]
            bas={}
            for x in bits:
                while x:
                    h=x.bit_length()-1
                    if h in bas: x^=bas[h]
                    else: bas[h]=x; break
            x=tb
            while x:
                h=x.bit_length()-1
                if h in bas: x^=bas[h]
                else: break
            res[Tset]=(x==0, len(V))
    return res
if __name__=="__main__":
    a=int(sys.argv[1]); t=int(sys.argv[2]); zeta=Fr(sys.argv[3]); b=int(sys.argv[4]); rho=int(sys.argv[5])
    r=run(a,t,zeta,b,rho)
    if r is None: print("abort"); sys.exit()
    bad=[T for T,(ok,_) in r.items() if not ok]
    print(a,t,zeta,b,rho,"criterion holds for all T:",not bad,"failures:",bad,flush=True)
