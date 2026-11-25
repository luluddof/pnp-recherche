# rho(D,y): Smith-weight tuple-count exponent / paper (D.15) exponent at one node, leading order.
# paper: min_th th*p + N*zeta2/th = 2 sqrt(zeta2 N p);  y := sqrt(zeta2 N/p).
# Smith: min_th th*p + N*sum_{k<D/2} f(th s_k), s_k = log2 D + 1 + 2k + 2 v2(k!), f(x)=-ln(1-e^-x).
import math
Z2 = math.pi**2/6
def v2f(k):  # v2(k!)
    s=0; 
    while k: k//=2; s+=k
    return s
def G(D, th, cache={}):
    key=D
    if key not in cache:
        c0 = int(math.log2(D)) + 1
        cache[key] = [c0 + 2*k + 2*v2f(k) for k in range(D//2)]
    return sum(-math.log1p(-math.exp(-th*s)) for s in cache[key] if th*s < 60)
def rho(D, y):
    # minimize th + (y^2/Z2) G(th) over th (golden search in log th)
    lo, hi = math.log(1e-7), math.log(50.0)
    g = lambda lt: math.exp(lt) + y*y/Z2*G(D, math.exp(lt))
    for _ in range(80):
        m1 = lo + (hi-lo)*0.382; m2 = lo + (hi-lo)*0.618
        if g(m1) < g(m2): hi = m2
        else: lo = m1
    return g((lo+hi)/2)/(2*y)
for la in (8, 11, 14, 17, 20):
    D = 2**la
    print(f"D=2^{la}: " + " ".join(f"y={y:g}:{rho(D,y):.4f}" for y in (1.0, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 1e-4)))
