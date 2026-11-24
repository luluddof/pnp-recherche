"""SAT search: is the crossing on N=4 (lambda=2) a GF(2)-sum of a row-rank-1 and b col-rank-1 functions?  (P_GF2(4) <= a+b)"""
import itertools, sys, time, threading
from pysat.solvers import Glucose4
def run(N,a,b,tl):
    nv=[0]
    def new():
        nv[0]+=1; return nv[0]
    cl=[]; terms=[]
    idx=list(itertools.product(range(N),repeat=2))
    for j in range(a+b):
        f={p:new() for p in idx}; g={p:new() for p in idx}; terms.append((j<a,f,g))
    for x,y,u,v in itertools.product(range(N),repeat=4):
        lits=[]
        for row,f,g in terms:
            p=new(); A=f[(x,y)] if row else f[(x,u)]; B=g[(u,v)] if row else g[(y,v)]
            cl += [[-p,A],[-p,B],[p,-A,-B]]; lits.append(p)
        rhs=int(x==v and y==u); k=len(lits)
        for signs in itertools.product((0,1),repeat=k):
            if sum(signs)%2!=rhs: cl.append([(-l if s else l) for l,s in zip(lits,signs)])   # forbid assignment 'signs'
    s=Glucose4(bootstrap_with=cl)
    t0=time.time()
    timer=threading.Timer(tl,lambda:s.interrupt()); timer.start()
    r=s.solve_limited(expect_interrupt=True); timer.cancel()
    out={True:"SAT",False:"UNSAT",None:"timeout"}[r]
    print("N=%d rows=%d cols=%d: %s (%.0fs)"%(N,a,b,out,time.time()-t0)); sys.stdout.flush()
    if r:
        m=set(l for l in s.get_model() if l>0)
        for row,f,g in terms: print("   ","row" if row else "col",[int(f[p] in m) for p in idx],[int(g[p] in m) for p in idx])
    s.delete(); return r
if __name__=="__main__":
    run(2,1,1,20); run(2,1,0,20); run(2,2,0,20)
    tl=int(sys.argv[1])
    for k in (4,5,6,7):
        for a in range((k+1)//2,k):
            run(4,a,k-a,tl)
