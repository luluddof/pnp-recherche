# Best achievable exponent/n of the (re-optimised) Claude recursion as a function of log2 m, density 4.26 and 1.
import math
rho=1/math.sqrt(2)
def hoeff(x,z):
    t=1.0
    for _ in range(100): t=max(1.0,(math.log(2*t)+2*x*math.log(2))/(2*z*z))
    return t
def levels(x0,z0,q,e0,ratio,L=20):
    e,x=e0,x0; out=[(0,x,e)]
    for l in range(L):
        z=z0*q**l; t=ratio(x,z)
        e=(rho+3*z)*e+(3/2**x if x<1000 else 0); x=2*(x+math.log2(t)); out.append((l+1,x,e))
    return out
configs={
 'paper (e0=1/2,x0=18,a=9)': [ (lambda: [(l, (3*18+6)*2**l-6*l-2*18-6, None) for l in range(20)]) ],
}
def best_e(logm, a, b, e0f, ratio):
    best=1e9
    for x0 in range(4,25):
        for k in range(2,10):
            for q in (0.5,0.6,0.7,0.8):
                for (l,x,e) in levels(x0,2.0**-k,q,e0f(x0),ratio):
                    if a*x+b<=logm: best=min(best,e)
    return best
def paper_e(logm):
    # paper recursion exactly
    e,x=0.5,18; best=1e9
    for l in range(20):
        if 9*x<=logm: best=min(best,e)
        z=2.0**(-l-5); e=(rho+3*z)*e+(3/2**x if x<1000 else 0); x=2*(x+3*l+18)
    return best
print(" log2 m |  paper: 2eW/n (r=4.26)  r=1   | optimised(C-E): r=4.26   r=1   | PPSZ 0.386, #SAT 0.716")
for lm in [20,40,64,100,200,400,1000,2000,5000,1e4,3e4,1e5,3e5,1e6]:
    ep=paper_e(lm); eo=best_e(lm,2,15,lambda x0:0.25+0.5/2**x0,hoeff)
    f=lambda e,w: (2*e*w if e<1e8 else float('nan'))
    print(f"{lm:>8g} | {f(ep,12.78):10.3f} {f(ep,3):10.3f} | {f(eo,11.78):10.3f} {f(eo,2):10.3f}")
