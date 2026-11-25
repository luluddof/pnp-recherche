# (A) check Smith(Per_k E) <= Smith(C_k E) termwise, faster Smith (vectorized pivot search). Usage: a t zeta k
import sys,itertools,time
from fractions import Fraction as Fr
import numpy as np
from lib12 import cauchy_block, toZ, v2, v2q, smith_vals
from levelm3 import perZ
NB=260;MOD=1<<NB
BIG=1<<(NB+5)
def smith_fast(A):
    M=np.array(A,dtype=object)%MOD; r,c=M.shape; out=[]
    for k in range(min(r,c)):
        sub=M[k:,k:]
        low=np.where(sub!=0, sub & (-sub), BIG)
        idx=np.argmin(low); ii,jj=divmod(int(idx),sub.shape[1])
        lv=int(low[ii,jj])
        if lv==BIG: out+=[NB]*(min(r,c)-k); break
        bv=lv.bit_length()-1; pi,pj=ii+k,jj+k
        M[[k,pi],:]=M[[pi,k],:]; M[:,[k,pj]]=M[:,[pj,k]]
        u=int(M[k,k])>>bv; M[k,:]=(M[k,:]*pow(u,-1,MOD))%MOD
        f=M[k+1:,k]>>bv
        M[k+1:,:]=(M[k+1:,:]-np.outer(f,M[k,:]))%MOD
        out.append(bv)
    return sorted(out)
if __name__=="__main__":
    a=int(sys.argv[1]); t=int(sys.argv[2]); zeta=Fr(sys.argv[3]); k=int(sys.argv[4])
    t0=time.time()
    E,c,T,B=cauchy_block(a,t,zeta,0); n=len(E)
    Ez=[[toZ(z,MOD) for z in r] for r in E]
    sE=smith_fast(Ez)
    subs=list(itertools.combinations(range(n),k))
    P=[[perZ([[Ez[r][u] for u in U] for r in R]) for U in subs] for R in subs]
    sp=smith_fast(P)
    sc=sorted(sum(x) for x in itertools.combinations(sE,k))
    ok=all(p<=q for p,q in zip(sp,sc))
    print(f"a={a} t={t} i={v2(t)} zeta={zeta} k={k} size={len(subs)}: (A) {ok}; min slack {min(q-p for p,q in zip(sp,sc))}; #ties {sum(p==q for p,q in zip(sp,sc))}; maxSmith {max(sp)}; {time.time()-t0:.0f}s",flush=True)
