# Counted App-B curve for the rank-m Heisenberg algebra (m colours): dims p_m(w) = [q^w] prod (1-q^n)^{-m}.
# R cost(s0) ~ P_m(2h+s0) (cumulative), C cost ~ P_m(2h-s0); N = p_m(h).  Uses log-space floats.
import math, sys
def logp_m(M, m):
    # exact big ints via recurrence w*p(w) = sum_{k=1}^w m*sigma(k) p(w-k)
    sig=[0]*(M+1)
    for d in range(1,M+1):
        for k in range(d,M+1,d): sig[k]+=d
    p=[1]+[0]*M
    for w in range(1,M+1):
        p[w]=sum(m*sig[k]*p[w-k] for k in range(1,w+1))//w
    return p
def run(h, m, out):
    p=logp_m(4*h+1, m); N=p[h]; L=math.log(N)
    P=[0]*(4*h+2); acc=0
    for w in range(4*h+1): acc+=p[w]; P[w]=acc
    rows=[]
    for s0 in range(2*h, -2*h-1, -1):
        th=math.log(P[2*h+s0])/L if 2*h+s0>=0 else 0
        T=math.log(P[2*h-s0])/L if 2*h-s0>=0 else 0
        rows.append((s0, th, T))
    # find jump: largest T among points with theta < 2 nearest the corner; min T with theta<=2
    ok=[(th,T) for s0,th,T in rows if th<=2]
    Tstar=min(T for th,T in ok) if ok else None
    bal=min(max(th,T) for s0,th,T in rows)
    end=math.log(P[4*h])/L
    print(f"h={h} m={m} kappa={m/h:.3f} log2N={L/math.log(2):.1f} theta_end=log P(4h)/log N={end:.3f}  min T s.t. theta<=2: {Tstar:.3f}  balanced={bal:.3f}", file=out, flush=True)
    return rows
if __name__=="__main__":
    out=sys.stdout
    for h in [60, 200]:
        for m in [1, 2, 4, max(1,h//20), h//5, h//2, h, 2*h]:
            run(h, m, out)
