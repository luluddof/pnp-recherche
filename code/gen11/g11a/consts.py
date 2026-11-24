import mpmath as mp
mp.mp.dps = 30
H = lambda x: -x*mp.log(x) - (1-x)*mp.log(1-x)
def kappa(b, e): return 3*(1-e)*H(b) - mp.log(2) - H(e)
for b, e in [(mp.mpf(1)/6, mp.mpf(1)/20), (mp.mpf(1)/8, mp.mpf(1)/50)]:
    k = kappa(b, e)
    print("beta", b, "eps", e, "H(b)", H(b), "H(e)", H(e), "3(1-e)H(b)", 3*(1-e)*H(b), "kappa", k, "ln2/kappa", mp.log(2)/k)
for den in [9, 60, 63]:
    rhs = 1/(den*mp.e**2)
    y = mp.findroot(lambda y: y*mp.log(mp.e/y) - rhs, 1e-4)
    print("den", den, "rhs", rhs, "eta", y)
c0 = (1-mp.e**-2)/4
c3 = c0*2/(mp.sqrt(2)*9)
print("c0", c0, "c3", c3, "c3 formula (1-e^-2)/(18 sqrt2)", (1-mp.e**-2)/(18*mp.sqrt(2)))
eta = mp.findroot(lambda y: y*mp.log(mp.e/y) - 1/(60*mp.e**2), 1e-4)
print("c3*eta*", c3*eta)
print("e^-6", mp.e**-6, "1/120", mp.mpf(1)/120, "14400 e^-6", 14400*mp.e**-6)
print("101/360>=1/6", mp.mpf(1)/3-mp.mpf(1)/20-mp.mpf(1)/360, "61/360", mp.mpf(2)/9-mp.mpf(1)/20-mp.mpf(1)/360)
# remark: L<=m/21, r<=m/600, N0<=m/600, u<=r, n0<=N0+3r
r = mp.mpf(1)/600; N0 = r; u = r; n0 = N0+3*r
print("remark: (r+u+n0)/3 =", (r+u+n0)/3, " mB>=", mp.mpf(2)/9-mp.mpf(1)/21-(r+u+n0)/3, " cross<=", mp.mpf(1)/21+r, " mA>=", mp.mpf(1)/3-mp.mpf(1)/21-(r+u+n0)/3)
print("e^-6 vs Delta/600 at Delta=2:", mp.e**-6, mp.mpf(2)/600)
# exact first moment (sharp trinomial) as independent check: sup over alpha, i, i' of
# H(alpha)/Delta... exponent per n:  ln-binomial(n,alpha n)/n + Delta * max_{(a,b,c)} [H3(a,b,c) + a ln alpha^3 + b ln (1-alpha)^3 + c ln(1-alpha^3-(1-alpha)^3)]
import numpy as np
def expo(Delta, grid=400):
    best = -1e9; arg = None
    for alpha in np.linspace(0.02, 0.98, 193):
        pa, pb = alpha**3, (1-alpha)**3; pc = 1-pa-pb
        Ha = -alpha*np.log(alpha)-(1-alpha)*np.log(1-alpha)
        # constraint a>=1/6, b>=1/6, c<=1/20, a+b+c=1
        for c in np.linspace(1e-6, 0.05, 26):
            a = np.linspace(1/6, 1-c-1/6, 200); b = 1-c-a
            val = -(a*np.log(a/pa) + b*np.log(b/pb) + c*np.log(c/pc))
            v = Ha + Delta*val.max()
            if v > best: best = v; arg = (alpha, c, a[val.argmax()])
    return best, arg
for D in [1.0, 1.2, 1.5, 1.766, 2.0, 3.0]:
    print("Delta", D, "sharp first-moment exponent/n (must be <0):", expo(D))
print("claimed crude exponent at Delta=2: ln2 - 2 kappa =", float(mp.log(2)-2*kappa(mp.mpf(1)/6, mp.mpf(1)/20)))
