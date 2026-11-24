from objs import *
from pred import comps
C1 = Copy([[2]], (Fr(1,2),)); C2=Copy([[2,0],[0,2]],(Fr(1,2),Fr(1,2)))
def F1(e,ha,hb):
    F=None
    for g in C1.S:
        t=kron(C1.B(g,ha,hb),C1.O(g,e)); F=t if F is None else madd(F,t)
    return F
def blockdiag_rank_tensor(E0,ha,hb):
    # rows: (ha-split, e-split) ; F2 = kron(F1(e1,ha1,hb1),F1(e2,ha2,hb2)) with index reorder (a1,m1',a2,m2') - rank unaffected by permutation
    tot=0
    for es in comps(E0,2):
        rowsplits=list(comps(ha,2)); colsplits=list(comps(hb,2))
        big=[]
        for hs in rowsplits:
            rowblk=None
            for ks in colsplits:
                blk=kron(F1(es[0],hs[0],ks[0]),F1(es[1],hs[1],ks[1]))
                rowblk=blk if rowblk is None else [r1+r2 for r1,r2 in zip(rowblk,blk)]
            big+=rowblk
        tot+=rankQ(big)
    return tot
for (e,ha,hb) in [(0,0,0),(2,0,0),(0,2,2),(1,1,1),(2,2,2)]:
    print(e,ha,hb,"tensor:",blockdiag_rank_tensor(e,ha,hb),"lattice:",flat_rank(C2,e,ha,hb))
