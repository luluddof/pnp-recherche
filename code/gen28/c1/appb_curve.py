# Counted App-B threshold curve with exact partition numbers: theta(s0) = log P_R(s0)/log N, T(s0) = log P_C(s0)/log N
# P_R(s0) = sum_{s<s0} sum_{i=0}^{2h} p(s+i) (row terms), P_C(s0) = sum_{s>=s0} sum_i p(2h-s-i); N = p(h).
import math, sys
def parts(M):
    p=[1]+[0]*M
    for n in range(1,M+1):
        k=1; s=0
        while True:
            g1=k*(3*k-1)//2
            if g1>n: break
            sg = 1 if k%2 else -1
            s+=sg*p[n-g1]
            g2=k*(3*k+1)//2
            if g2<=n: s+=sg*p[n-g2]
            k+=1
        p[n]=s
    return p
def lg(x): return math.log(x) if x>0 else float('-inf')
for h in [int(a) for a in sys.argv[1:]]:
    p=parts(4*h+2); N=p[h]; L=math.log(N)
    P=lambda w: p[w] if 0<=w<=4*h else 0
    # cumulative: R(s0): sum_{s=-2h}^{s0-1} sum_{i=0}^{2h} p(s+i)
    rowc={s: sum(P(s+i) for i in range(2*h+1)) for s in range(-2*h,2*h+1)}
    colc={s: sum(P(2*h-s-i) for i in range(2*h+1)) for s in range(-2*h,2*h+1)}
    R=0; Rs={}
    for s0 in range(-2*h,2*h+2):
        Rs[s0]=R
        if s0<=2*h: R+=rowc[s0]
    C=0; Cs={}
    for s0 in range(2*h+1,-2*h-1,-1):
        if s0<=2*h: C+=colc[s0]
        Cs[s0]=C
    Cs[2*h+1]=0
    print(f"h={h} N=p(h)={N:.4g}  log2 N={L/math.log(2):.1f}")
    print("  w=2h-s0+... : s0, eta=T, theta(capped at 2), 2-theta, eta^2/4, circle-theta")
    for s0 in sorted(set([2*h+1-w for w in [1,2,3,5,8,12,20,30,50,80,120,200,400,800,1600] if 2*h+1-w>-2*h]+[0])) [::-1]:
        T=lg(Cs[s0])/L; th=min(2,lg(Rs[s0])/L) if Rs[s0]>0 else 0
        sig=(s0-0.5)/h
        circ=math.sqrt(max(0,2+sig))
        print(f"  s0={s0:6d} eta={T:.4f} theta={th:.4f} delta={2-th:.5f} eta^2/4={T*T/4:.5f} [asympt sqrt(2+s/h)={circ:.4f}, T_asym={math.sqrt(max(0,2-sig)):.4f}]")
