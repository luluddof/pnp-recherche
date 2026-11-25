# Flexible explicit model: rows B.. of the ablation.  See appd_explicit.py header for conventions.
import math
from appd_explicit import ladd, PIB, LOGN_OFF, W_INC, W_PAPER, lam_lin, eps_legal, threshold
from scipy.optimize import minimize_scalar, minimize
def run(lam, us, cs=None, w=W_PAPER, kappa='24', pack='paper', Dmode='s09', amode='paper', r=64,
        leaf='paper', smith=False, Delta0=64, S2root=17.0, verbose=False, leaffac=1.0):
    """us: list of u_d (depth schedule); recursion stops at leaf rule.
       kappa: '24' (8hD,h=3M/u), '18' (6hD,h=3M/u), 'lin' (exact loss of linear taper, h=M/(u-1)), or float
       pack: 'paper' (nu=6, child 3M/32 under (D.6)) | 'slack' (nu=2P(1-1/r), child 2PM/r, P from Lemma D.3 proof)
       Dmode: 's09' (Delta = s^0.9, S2 <= 2 M Delta) | 'adapt' (S2 <= M(Delta0 + 2 L_parent))
       amode: 'paper' (a_d = a_{d-1}+2u_d) | 'slack' (a_d = a_{d-1}+u_d+c_d)
       leaf: 'paper' (M <= s/floor(sqrt lam)) | 'depth' (exactly len(us) levels)"""
    logn=lam-LOGN_OFF; LM=math.log2(w); Lp=0.0; a_prev=0; a_pp=0; E=0.0; info=[]
    Hn = 2**LOGN_OFF/math.floor(math.sqrt(lam))
    d=0
    while True:
        if leaf=='paper' and 2**LM <= Hn: break
        if leaf=='depth' and d==len(us): break
        if d>=len(us): return float('inf')
        u=us[d]
        if u<2: return float('inf')
        a = a_prev+2*u if amode=='paper' else a_prev+u+cs[d]
        logM=logn+LM
        if pack=='paper':
            Dlog=2*a*math.log2(r)
            # (D.6) check at smallest active M (H+1) -- (paper) ; here at this M
            if 1+2*Dlog+0.9*lam > logM+1e-9 or a+Dlog>logM: return float('inf')
            nufac=6.0; shrink=6/r
        else:
            Dlog=2*a*math.log2(r)
            if Dmode=='s09': S2o=1+0.9*lam
            else: S2o = math.log2(S2root) if d==0 else math.log2(Delta0+2**(a_prev+1)+(2**(a_pp+1) if d>=2 else 0))
            t1=Dlog+0.5*S2o-0.5*logM; t2=a+Dlog-logM
            if t1>40 or t2>40: return float('inf')
            P=1+2**t1+2**t2
            nufac=2*P*(1-1/r); shrink=2*P/r
        if smith: nufac/=4
        best=0.0
        for i in sorted(set([0,1,2,u-2,u-1])):
            if i<0 or i>=u: continue
            if kappa=='24': Lam=24/u; ev=8*(a+1)
            elif kappa=='18': Lam=18/u; ev=4*(a+1)
            elif kappa=='lin': Lam=lam_lin(u,a,i,1/(u-1)); ev=(u-i)/(u-1)*(a+1-i)
            else: Lam=float(kappa)/u; ev=(u-i)/(u-1)*(a+1-i)
            Dl=a-i
            best=max(best, Lam + 2**(Lp-LM-Dl) + ev*2**(-Dl))
        E += PIB*math.sqrt(2/3*nufac*best)*2**LM
        info.append((u,a,round(nufac,3),round(2**LM,3)))
        Lam0 = {'24':24/u,'18':18/u}.get(kappa, None)
        if Lam0 is None: Lam0 = lam_lin(u,a,0,1/(u-1)) if kappa=='lin' else float(kappa)/u
        Lp=ladd(Lp, math.log2(Lam0)+LM+a); Lp=ladd(Lp, math.log2(8*(a+1))+LM)
        LM+=math.log2(shrink); a_pp=a_prev; a_prev=a; d+=1
    E += (Hn if leaf=='paper' else leaffac*2**LM)
    if verbose: print("    levels (u,a,nufac,M/n):",info,"leaf/n=%.3g"%(Hn if leaf=='paper' else 2**LM))
    return E
def paper_sched(lam, eps, depth=12):
    return [math.floor(eps*lam/2**d) for d in range(depth)]
