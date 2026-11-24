# G14f exp2: leading-order exponents (in units of N = n/2, log2) of every encoding on the bipartite benchmark,
# at the benchmark scale e = 1/d  (d ~ sqrt(log n)).  Constants of the paper's O(.) set to the explicit ones where
# they are explicit (U formulas), cost of one network = e*U.
from math import log2, e as E, ceil, exp, lgamma, log

def H(p):
    return 0.0 if p <= 0 or p >= 1 else -p*log2(p) - (1-p)*log2(1-p)

print("A. whole-sum encodings: exponent / N  (trivial enumeration of one side = 1)")
print("   d | clause-net 2e(d-1) | 3-CNF chain | Fourier Z_{d+1} | Z_2^d characters | IS 2-CNF direct | Prop.R lead. 2e*log2(2d) (x const>=28) | target 1/d")
for d in [4, 8, 16, 32, 64, 256, 1024]:
    e = 1.0/d
    clause = e*2*(d-1)
    # chain: prefix-OR gadget: y_i <-> y_{i-1} v x_i ; aux vars d-1 per factor, each aux occurs in 2 gates (+1 weight factor),
    # tests = sum_var (occ-1): originals N(d-1), aux N(d-1)*1  -> U = 2*2N(d-1)
    chain = e*2*2*(d-1)
    four = e*2*(d-1)*ceil(log2(d+1))
    char = e*2*d            # one bit Gamma_v per edge, shared by 2 clusters
    direct = e*2*d          # U=Theta(Nd) (G8b): 1 bit per crossing edge, both ends
    propR = e*2*log2(2*d)
    print(f"{d:5d} | {clause:8.3f} | {chain:8.3f} | {four:8.3f} | {char:8.3f} | {direct:8.3f} | {propR:8.4f} | {1/d:8.4f}")

print()
print("B. class |T| = tau N/d (T in R sparse, A subset of L\\N(T) free): exponent / (N/d)")
print("   d  tau | enumerate T | T broadcast (Lemma BC) e*2tau*log(3d/tau) | aggregation-model LB e*(d-2)H(tau/d)*d | cut-rank LB e*d*H(tau/d) | ratio best/target(=1)")
for d in [16, 64, 256, 1024, 2**16]:
    e = 1.0/d
    for tau in [0.25, 1.0, 2.0]:
        enum = d*H(tau/d)
        bc = d*e*2*tau*log2(3*d/tau)
        agg = d*e*(d-2)*H(tau/d)
        cut = d*e*H(tau/d)*d/d
        print(f"{d:6d} {tau:4.2f} | {enum:8.3f} | {bc:8.3f} | {agg:8.3f} | {cut:8.4f} | {min(enum,bc):8.3f}")

print()
print("C. refined classes (|T|,|A|) = (tau N/d, alpha N/d): min of {T broadcasts, A broadcasts, Thm B (const 1), enumerate both}, e=1/d;")
print("   maximised over feasible (tau, alpha) [alpha <= d e^{-tau} OR tau <= d e^{-alpha}]; exponent/(N/d)")
def cost(d, tau, al):
    e = 1.0/d
    bT = e*2*tau*log2(3*d/tau)*d
    bA = e*2*al*log2(3*d/al)*d if al < 3*d else 1e9
    h = d*(H(tau/d) + H(min(al/d, 0.5)))          # log|class| / (N/d)
    thB = e*h*log2(4*d*d/max(h,1e-9))               # e h log(4m/h), m = N d, h in units N/d
    en = h
    return min(bT, bA, thB, en), (bT, bA, thB, en)
for d in [16, 64, 256, 1024, 2**16]:
    best = (0, None)
    grid = [2**(i/4) for i in range(-40, 80) if 2**(i/4) <= d/2]
    for tau in grid:
        for al in grid:
            if al > d*exp(-tau) and tau > d*exp(-al): continue   # feasible if EITHER side is typical
            c, parts = cost(d, tau, al)
            if c > best[0]: best = (c, (tau, al, parts))
    c, (tau, al, parts) = best
    print(f"   d={d:6d}: max-min cost = {c:7.3f} (N/d)  at tau={tau:.3g}, alpha={al:.3g} (alpha/d={al/d:.3g});  log2 d = {log2(d):.1f};"
          f" parts bT,bA,thB,enum = " + ", ".join(f"{x:.2f}" for x in parts))

print()
print("D. balanced-cut log-rank of F(A)=2^{|R\\N(A)|} (upper estimate via traces), per N:  max_tau min{H(tau/d), 0.5*H(e^-tau)}")
for d in [16, 64, 256, 1024, 2**16]:
    best = max(min(H(t/d), 0.5*H(exp(-t))) for t in [0.01*i for i in range(1, 4000)])
    print(f"   d={d:6d}: <= {best:.4f} N   ( (log2 d)^2/d = {log2(d)**2/d:.4f} ; single-cluster cut ranks sum to N )")
