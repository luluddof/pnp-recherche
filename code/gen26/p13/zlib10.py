# G26-P10: own transcription of the deviation terms (from the lemma statements, not from P5's code).
import math, numpy as np
LN2 = math.log(2)
def _sumsqrt(t):            # sum_{j=0}^{t-1} sqrt(j): exact for t <= 2^20, else integral upper bound (2/3) t^1.5
    t = float(t)
    if t <= 2**20:
        return float(np.sqrt(np.arange(int(t))).sum())
    return (2/3)*t**1.5
def zF(lt, x, d1=0.33, d2=0.33, d3=0.33, mcd=True):
    t = 2.0**lt; K = 2.0**min(x, 60)
    L = math.log(4/d1) + (lt + 2*x)*LN2; ell = math.log(2/d3) + (x + lt)*LN2
    L2 = math.log(K*(K-1)/(2*d2)) if x < 60 else 2*x*LN2 - math.log(2*d2)
    s2 = math.sqrt(2*(1 + 1/(K-1))) if x < 60 else math.sqrt(2)
    m = min(math.sqrt(L/2), 0.5 + math.sqrt(ell/2)) if mcd else math.sqrt(L/2)
    # V' = sum_{k=1}^t (s2 (m sqrt(k-1) + 1/2) + 1/2)^2 = s2^2 m^2 t(t-1)/2 + s2 m (s2+1) sum sqrt(k-1) + t ((s2+1)/2)^2
    V = s2*s2*m*m*t*(t-1)/2 + s2*m*(s2+1)*_sumsqrt(t) + t*((s2+1)/2)**2
    b = 4*math.sqrt((t-1)*L/2) + 3
    s = b*L2/3 + math.sqrt((b*L2/3)**2 + 2*V*L2)
    return 2/t + s/t**2 + (1/(K-1) if x < 60 else 0.0)
def zR2(lt, x, d1=0.5, d2=0.49):
    t = 2.0**lt; K = 2.0**min(x, 60)
    L = math.log(2/d1) + 2*(lt + x)*LN2; L2 = 2*x*LN2 - math.log(2*d2)
    # Sig = sum_{k<t} (4 sqrt(k L/2) + 1)^2 = 8L t(t-1)/2 + 8 sqrt(L/2) sum sqrt(k) + t
    Sig = 4*L*t*(t-1) + 8*math.sqrt(L/2)*_sumsqrt(t) + t
    return 4/t + math.sqrt(2*Sig*L2)/t**2 + (1/(K-1) if x < 60 else 0.0)
def sdS2(lt, x):            # K^2 sd(Y)/N  (Lemma S2 bound)
    t = 2.0**lt; K = 2.0**min(x, 60)
    s2 = math.sqrt(2*(1 + 1/(K-1))) if x < 60 else math.sqrt(2)
    # sum_{k=1}^t (s2 (sqrt(k-1)/2 + 1/2) + 1/2)^2
    V = s2*s2*t*(t-1)/8 + s2*(s2+1)/2*_sumsqrt(t) + t*((s2+1)/2)**2
    return (1/math.sqrt(2) + math.sqrt(V)/t)/t
def zS2(lt, x, factor=1.0):  # factor 1 = corrected (both diagonals padded), 0.5 = P5's
    K = 2.0**min(x, 60)
    return factor*sdS2(lt, x) + (1/(K-1) if x < 60 else 0.0)
def r3terms(lg, version='R3', dA=0.33, dG=0.33, dGam=0.33):
    """exact discrete sums (no integral bounds) for g <= 2^24, else fall back to fine trapezoid on the exact summand.
    Gamma: c_i = 4 G_i / w_{i+1}, G_i = R sqrt(i LG/2) [x sqrt(1-(i-1)/(g-1)) Serfling], S_j = sum_{i<j} c_i^2,
    dev_j = w_j sqrt(2 S_j Lam) (Azuma) or w_j sqrt(S_j Lam/2) (interval Hoeffding-Azuma). Returns (cA, Gam) per N."""
    g = 2.0**lg
    ser = version in ('R3s', 'R3pp'); half = version in ('R3h', 'R3pp')
    cA = math.log(g/dA)*(g-1)/(4*(g-2)**2)*((0.5 + 1/g) if ser else 1.0)
    R = 1/g + 1/(g-2); LG = math.log(2*g*g/dG); Lam = math.log(g/dGam)
    n = int(g//2) + 1
    if g <= 2**24:
        i = np.arange(0, n, dtype=float)
    else:
        i = np.linspace(0, g/2, 2_000_001)
    w1 = (g - i - 1)*(g - i - 2)/(g*(g - 1))
    G2 = R*R*i*LG/2*(np.maximum(1 - (i - 1)/(g - 1), 0) if ser else 1.0)
    c2 = 16*G2/w1**2
    if g <= 2**24:
        S = np.concatenate([[0.0], np.cumsum(c2)])[:-1]     # S_j = sum_{i<j}
        j = i
    else:
        di = np.diff(i); S = np.concatenate([[0.0], np.cumsum((c2[1:] + c2[:-1])/2*di)])*1.0005; j = i
    w = (g - j)*(g - j - 1)/(g*(g - 1))
    dev = w*np.sqrt((S*Lam/2) if half else (2*S*Lam))
    return cA, dev.max()/g
