from objs import *
C1 = Copy([[2]], (Fr(1,2),))
for (e,h) in [(0,0),(2,2),(4,4),(6,6),(4,0),(0,4),(6,2),(2,6)]:
    tot,bl=flat_rank(C1,e,h,h,1,True)
    print(e,h,"rkF",tot,"blocks(rows,cols,rank)",sorted(bl,reverse=True),flush=True)
print("dims",[len(C1.basis(k)) for k in range(10)])
