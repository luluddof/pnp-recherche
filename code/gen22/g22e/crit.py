import math,sys
from chars import *
SQ2=math.sqrt(2)
def cum(d): 
    s=0;D=[]
    for x in d: s+=x; D.append(s)
    return D
def L(D,E,den=1):  # log dim M_{<=E}; E real, den = grading denominator
    i=int(math.floor(E*den+1e-9)); i=min(i,len(D)-1); return math.log(D[i])
def report(name,d,den=1,c=None,ceff=None,ground=1.0):
    D=cum(d); Emax=(len(D)-1)/den
    print("==",name,"Emax",Emax)
    # (a) Cardy check
    if ceff is not None:
        A=2*math.pi*math.sqrt(ceff/6)
        for E in [Emax/16,Emax/4,Emax]:
            print("  E=%7.1f  L=%9.3f  A*sqrt(E)=%9.3f  ratio=%.4f" % (E,L(D,E,den),A*math.sqrt(E),L(D,E,den)/(A*math.sqrt(E))))
    # (b) criterion (*) for M=V, without log terms: R(w)=2L(2w)/min(L(w),L(4w)/2)
    worst=(9,None); g=[]
    w=1.0/den
    while 4*w<=Emax:
        Lw,L2,L4=L(D,w,den),L(D,2*w,den),L(D,4*w,den)
        if Lw>0:
            R=2*L2/min(Lw,L4/2)
            if R<worst[0]: worst=(R,w)
        w+=1.0/den
    print("  min_w R(w) = %.4f at w=%s  (barrier 2sqrt2=%.4f)"%(worst[0],worst[1],2*SQ2))
    # g(E)=L/sqrt(E) monotone?
    gs=[(E/den,math.log(D[E])/math.sqrt(E/den)) for E in range(1,len(D))]
    ch=sum(1 for i in range(2,len(gs)) if (gs[i][1]-gs[i-1][1])*(gs[i-1][1]-gs[i-2][1])<0)
    print("  L/sqrtE: first %.3f, last %.3f, sign changes of increments: %d"%(gs[0][1],gs[-1][1],ch))
    return worst
if __name__=="__main__":
    N=int(sys.argv[1]) if len(sys.argv)>1 else 400
    report("Heisenberg (partitions)",heis(N),1,1,1)
    report("Heisenberg rank 8",heis(N,8),1,8,8)
    report("NS fermion n=1 (half units)",NSferm(2*N),2,0.5,0.5)
    report("NS fermions n=8",NSferm(2*N,8),2,4,4)
    report("E8 lattice",E8(N),1,8,8)
    report("Leech lattice VOA",leech(min(N,300)),1,24,24)
    report("Moonshine V#",J_dims(min(N,300)),1,24,24)
    report("Yang-Lee vacuum M(2,5)",RR(N,{2,3}),1,-22/5,2/5)
    report("Yang-Lee module h=-1/5",RR(N,{1,4}),1,-22/5,2/5)
    report("M(2,7) vacuum",RR(N,{3,4},7),1,-68/7,4/7)
    report("symplectic fermions (prod(1+q^k)^2)",symplf(N),1,-2,1)
    report("Ramond n=8 (no ground factor)",ramond(N,8),1,4,4)
