import numpy as np, itertools
def allx(M): return np.array(list(itertools.product((-1,0,1),repeat=M)),dtype=np.int64)
def Aexact(wires,n,parts):
    """wires: list (h,t) vertices; parts: list of arrays lab[v] in 0..r-1. returns A'' = max_x min_i b_i and per-x b's"""
    M=len(wires); D=np.zeros((M,n),dtype=np.int64)
    for e,(h,t) in enumerate(wires): D[e,h]+=1; D[e,t]-=1
    X=allx(M); Y=X@D
    bs=[]
    for lab in parts:
        r=lab.max()+1; P=np.zeros((n,r),dtype=np.int64); P[np.arange(n),lab]=1
        bs.append(np.abs(Y@P).sum(1)//2)   # integer: sum_g |y(g)| is even
    B=np.stack(bs,1)
    return int(B.min(1).max()), B, X
def Mi(wires,parts): return [sum(lab[h]!=lab[t] for h,t in wires) for lab in parts]
