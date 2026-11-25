# Explicit operation-count model of OpenAI App. D for #k-SAT at clause density alpha, s = n+m, lambda = log2 s.
# Exponent (bits, divided by n) of the worst root-to-leaf path:  E/n = sum_d f_d/n + M_leaf/n.
# f_d = (pi/ln2) * sqrt(2/3 * n_t p_t) at the worst node (the (D.15) bound 2 sqrt(zeta(2) n p)).
# Polynomial factors s^{O(1)} and the packet-entry cost exp(O(L log s)) are ignored (2^{O(lambda)} << 2^n here).
import math, sys
import numpy as np
from scipy.optimize import minimize
LN2=math.log(2); PIB=math.pi/LN2
def ladd(x,y):  # log2(2^x+2^y)
    m=max(x,y); return m+math.log2(2**(x-m)+2**(y-m))
K=3; ALPHA=4.26
LOGN_OFF=math.log2(1+ALPHA)          # log2 n = lambda - LOGN_OFF
W_INC=K*ALPHA                        # incidence network: M0 = k m  (12.78 n)
W_PAPER=2*(K*ALPHA+ (K*ALPHA-3))     # trees+subdivision: 2(km + sum(d_v-3)^+) ~ 45.12 n (Poisson(12.78), P(d<3)~0)
def lam_lin(u,a,i,hfac):   # Lambda_i / M for q(i)=h(u-i), h = hfac*M/u  (paper hfac=3; minimal: h=M/(u-1))
    h = hfac
    return h*sum(j/2**j*min(j,u-i) for j in range(1,min(a-i,200)+1))

# ---------------- paper as written (and simple variants) ----------------
def paper(lam, eps=2**-17, w=W_PAPER, kappa_mode='paper', nu=6.0, shrink=3/32, sched='paper', r=64, verbose=False):
    logn = lam-LOGN_OFF
    M=[w]               # M_d / n along worst path
    H = 2**(LOGN_OFF)/math.floor(math.sqrt(lam))   # H/n = s/(n floor sqrt lam)
    d=0; E=0.0; a_prev=0; us=[]
    while M[-1] > H:
        u = math.floor(eps*lam/2**d)
        if u < 2: return float('inf'), None
        a = a_prev + 2*u
        if kappa_mode=='paper':   X = 24/u                              # p_t ~ 8 h D, h = 3M/u
        elif kappa_mode=='loss6': X = 6*3/u                          # p_t ~ 6hD (true loss bound), h=3M/u
        elif kappa_mode=='hmin':  X = lam_lin(u,a,0,1/(u-1))         # h = M/(u-1), budget = exact loss
        E += PIB*math.sqrt(2/3*nu*X)*M[-1]
        us.append(u); a_prev=a
        M.append(M[-1]*shrink); d+=1
    E += H          # leaves (M <= H) by brute force
    if verbose: print("  depths",d,"u_d",us,"leaf H/n=%.3g"%H)
    return E, us
def eps_legal(lam, w=W_PAPER, r=64, delta_exp=0.9, shrink=3/32):
    # largest eps s.t. (D.6) holds at every possible active call (M >= H+1) for the paper's schedule
    s_log=lam; H_log = lam - math.log2(math.floor(math.sqrt(lam)))
    logn=lam-LOGN_OFF; M0_log = logn+math.log2(w)
    dmax=0
    while M0_log + (dmax+1)*math.log2(shrink) > H_log: dmax+=1
    def ok(eps):
        a=0
        for d in range(dmax+1):
            u=math.floor(eps*lam/2**d); a+=2*u
        Dlog=2*a*math.log2(r)
        return (1+2*Dlog+delta_exp*s_log + 1e-9 <= H_log) and (a+Dlog<=H_log) and (a+1<=delta_exp*s_log)
    lo,hi=0.0,1.0
    for _ in range(60):
        mid=(lo+hi)/2
        if ok(mid): lo=mid
        else: hi=mid
    return lo, dmax

# ---------------- improved explicit algorithm ----------------
def improved(lam, x, D, r=64, w=W_INC, smith=False, kappa=None, Delta0=64, S2root=17.0, verbose=False, lp=False):
    # x = [u_0..u_{D-1}, c_0..c_{D-1}]  (u real -> floor, c real -> ceil >=0)
    logn = lam-LOGN_OFF
    us=[max(2,math.floor(v)) for v in x[:D]]; cs=[max(0,math.ceil(v)) for v in x[D:2*D]]
    LM=math.log2(w); Lp=0.0    # log2(M_d/n), log2(p_inh/n)  (p_root = n+1)
    a_prev=0; E=0.0; info=[]
    for d in range(D):
        u=us[d]; a=a_prev+u+cs[d]
        logM=logn+LM
        Dlog=2*a*math.log2(r)
        S2over = math.log2(S2root) if d==0 else math.log2(Delta0+2**(a_prev+1))
        t1=Dlog+0.5*S2over-0.5*logM; t2=a+Dlog-logM
        if t1>60 or t2>60: return float('inf'), None
        P=1+2**t1+2**t2
        nufac=2*P*(1-1/r)
        if smith: nufac/=4
        best=0.0
        for i in ([0,u-1] if u>60 else range(u)):
            if kappa is not None: Lam = kappa/u if i==0 else kappa/u
            else: Lam=lam_lin(u,a,i,1/(u-1))
            if lp: Lam = Lam*min(1.0, 1.0)  # LP gain <0.3% at u>=300: not modelled
            q=(u-i)/(u-1)    # q_t/M
            ev = q*(a+1-i)                       # evaluation loss / M
            Dl = a-i
            term = Lam + 2**(Lp-LM-Dl) + ev*2**(-Dl)
            best=max(best,term)
        E += PIB*math.sqrt(2/3*nufac*best)*2**LM
        info.append((u,a,cs[d],round(P,4),round(2**LM,4)))
        # inherited precision for children: max_t p_t (i=0, D=L=2^a)
        Lam0 = (kappa/u) if kappa is not None else lam_lin(u,a,0,1/(u-1))
        Lp = ladd(Lp, math.log2(Lam0)+LM+a); Lp=ladd(Lp, math.log2((u/(u-1))*(a+1))+LM)
        LM = LM + math.log2(2*P/r)
        a_prev=a
    E += 2**LM   # leaves: brute force on M_D edges
    if verbose: print("   levels (u,a,c,P,M/n):",info," leaf M/n=%.3g"%2**LM)
    return E, info

def opt_improved(lam, D, r=64, **kw):
    # schedule u_d ∝ M_d^{2/3}: parametrise by u0, ratio, c's
    def f(y):
        u0=math.exp(y[0]); rho=math.exp(y[1]); c=y[2:]
        x=[u0*rho**d for d in range(D)]+list(c)
        return improved(lam,x,D,r=r,**kw)[0]
    best=None
    for g0 in (lam/30, lam/60, lam/120):
        for gr in (0.1,0.2,0.3):
            y0=[math.log(g0),math.log(gr)]+[12.0]*D
            res=minimize(f,y0,method='Nelder-Mead',options={'maxiter':3000,'xatol':1e-4,'fatol':1e-9})
            if best is None or res.fun<best.fun: best=res
    y=best.x; u0=math.exp(y[0]); rho=math.exp(y[1])
    x=[u0*rho**d for d in range(D)]+list(y[2:])
    return best.fun, x

def threshold(Efun, tau, lo=1e2, hi=1e14):
    # smallest lambda with E(lam)/n < tau (E/n decreasing in lam), geometric bisection
    for _ in range(60):
        mid=math.sqrt(lo*hi)
        if Efun(mid) < tau: hi=mid
        else: lo=mid
        if hi/lo<1.002: break
    return hi
