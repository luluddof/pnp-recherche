# Twisted sectors (orbifold): M = g-twisted V_L-module (fixed-point-free g, conformal grading => S = {0}, only the vacuum channel).
# N <= min(2 D_V^conf(H/2), D_M(E0)) is then RIGOROUS (dichotomy branch (i) is the only branch).  Per-copy criterion.
from scan import *
def phi_tw(ground, T, dims, r):
    out=[]
    for b in bs:
        s=math.log(ground)
        for j,d in dims.items():
            k=np.arange(0,4000); s-= d*np.log1p(-np.exp(-b*(k+j/T))).sum()
        out.append(s)
    return np.array(out)
A1=[[2]]; A2=[[2,-1],[-1,2]]; D4=[[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]]
E8=[[2,-1,0,0,0,0,0,0],[-1,2,-1,0,0,0,0,0],[0,-1,2,-1,0,0,0,-1],[0,0,-1,2,-1,0,0,0],[0,0,0,-1,2,-1,0,0],[0,0,0,0,-1,2,-1,0],[0,0,0,0,0,-1,2,0],[0,0,-1,0,0,0,0,2]]
cases=[("A1, g=-1, ground 1",A1,1,2,{1:1},60),("A2, g=-1, ground 2",A2,2,2,{1:2},40),("A2, g=order 3, ground 1",A2,1,3,{1:1,2:1},40),
       ("D4, g=-1, ground 4",D4,4,2,{1:4},20),("D4, g=order 3 fpf, ground 3",D4,3,3,{1:2,2:2},20)]
for name,G,gr,T,dims,R2 in cases:
    r=len(G); PM=phi_tw(gr,T,dims,r); PV=phi(G,[0]*r,[0]*r,R2)
    grid=np.exp(np.linspace(math.log(0.03*r),math.log(3.0*r),38))
    U=lambda u,v: rate(PM,2*u+2*v)/2
    c=crit(PM,[lambda u,v:rate(PM,u),lambda u,v:rate(PV,v),U],grid)
    own=crit(PM,[lambda u,v:rate(PM,u),lambda u,v:rate(PM,v),U],grid)
    print("%-32s TRUE (pairs in V conf) %.4f at e=%.2f v=%.2f | M as its own anchor (not a valid bound) %.4f at e=%.2f v=%.2f"%(name,*c,*own),flush=True)
# E8: theta via E4 is cheaper: phi_V(b) = log E4(e^-b) - 8 sum log(1-e^-bk)
k=np.arange(1,4000)
def phiE8():
    out=[]
    for b in bs:
        n=np.arange(1,3000); s3=np.array([sum(d**3 for d in range(1,m+1) if m%d==0) for m in range(1,200)])
        q=math.exp(-b); e4=1+240*sum(s3[m-1]*q**m for m in range(1,200)) if b>0.3 else None
        if e4 is None:  # modular: E4(q) = (2pi/b)^4 E4(q~), q~=exp(-4pi^2/b)
            qt=math.exp(-4*math.pi**2/b); e4=(2*math.pi/b)**4*(1+240*qt)
        out.append(math.log(e4)-8*np.log1p(-np.exp(-b*k)).sum())
    return np.array(out)
PV=phiE8()
for name,gr,T,dims in [("E8, g=-1, ground 16",16,2,{1:8}),("E8, g=order 3 fpf, ground 9",9,3,{1:4,2:4})]:
    PM=phi_tw(gr,T,dims,8); grid=np.exp(np.linspace(math.log(0.2),math.log(24),38))
    U=lambda u,v: rate(PM,2*u+2*v)/2
    c=crit(PM,[lambda u,v:rate(PM,u),lambda u,v:rate(PV,v),U],grid)
    own=crit(PM,[lambda u,v:rate(PM,u),lambda u,v:rate(PM,v),U],grid)
    print("%-32s TRUE (pairs in V conf) %.4f at e=%.2f v=%.2f | M as its own anchor %.4f at e=%.2f v=%.2f"%(name,*c,*own),flush=True)
