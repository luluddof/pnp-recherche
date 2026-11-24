import numpy as np, networkx as nx
from anat import *
def blocks(N,R):
    g=nx.Graph()
    for a in range(N*N):
        for b in range(N*N):
            if R[a,b]: g.add_edge(('r',a),('c',b))
    out=[]
    for comp in nx.connected_components(g):
        rows=sorted(a for t,a in comp if t=='r'); cols=sorted(b for t,b in comp if t=='c')
        sub=R[np.ix_(rows,cols)]
        out.append((rows,cols,sub))
    return out
if __name__=="__main__":
  for f in files:
    N,G,Gp=load(f); R,C=Rmat(N,G,Gp)
    print("=====",f.split('/')[-1])
    for rows,cols,sub in sorted(blocks(N,R),key=lambda t:-len(t[0])):
        print(" rows(x,y):",[(a//N,a%N) for a in rows]," cols(u,v):",[(b//N,b%N) for b in cols]," rank",np.linalg.matrix_rank(sub), "full" if (sub!=0).all() else "sparse")
        print("   ",sub.tolist())
