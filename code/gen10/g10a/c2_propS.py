import math, random
import numpy as np
ln = math.log
# (1) constant in (iii): rho = g*r <= 10 sqrt(b ln n) n/sqrt(s) ?
def iii(n, b, s):
    g = math.ceil(n/s); r = math.ceil(9*math.sqrt(b*s*ln(n)))
    return g*r, 10*math.sqrt(b*ln(n))*n/math.sqrt(s), g, r
bad = []
for n in (10**6, 10**9, 10**12):
    for frac in (n**(-1/3), 0.05, 0.1, 0.12, 0.2, 0.34, 0.51, 0.6, 0.9, 1.0):
        s = frac*n
        rho, bd, g, r = iii(n, 3, s)
        if rho > bd: bad.append((n, round(frac,4), g, r, rho, round(bd)))
print("(iii) violations (n, s/n, g, r, rho, bound):")
for x in bad: print("  ", x)
# (2) r <= 6s needs s >= ~2.25 b ln n; hypothesis is s >= 15 ln n
for b in (3, 6, 10, 20):
    n = 10**9; s = 15*ln(n)
    r = math.ceil(9*math.sqrt(b*s*ln(n)))
    print(f"b={b}: s=15 ln n={s:.0f}, r={r}, 6s={6*s:.0f}, r<=6s: {r<=6*s}, r>2s (trivial case): {r>2*s}")
# (3) union bound arithmetic
n=1000; b=3
X = math.ceil(n**3*ln(b)) + n**2
print("log P(fail, Hoeffding+union) <=", n*ln(b) - 2*X/n**2, "; R =", math.ceil(2*X/n), " q+1 <=", math.ceil((3*n*math.ceil(2*X/n)+1)/X))
# (4) simulate the construction at one point: random omega vs single-point adversarial omega
rng = np.random.default_rng(0)
def build(n):
    # random 3-regular formula: n clauses, each var 3 times (configuration model, retry on repeats)
    while True:
        p = rng.permutation(np.repeat(np.arange(n), 3)).reshape(n, 3)
        if all(len(set(e)) == 3 for e in p): return p
n = 3000; b = 3; s = round(n**(2/3)); g = math.ceil(n/s)
cl = build(n)
occ = [[] for _ in range(n)]
for e, vs in enumerate(cl):
    for v in vs: occ[v].append(e)
lab = rng.integers(0, g, n); ue = cl[np.arange(n), rng.integers(0, 3, n)]; alpha = lab[ue]
a = np.bincount(lab, minlength=g)
mx = 0
for _ in range(300):
    om = np.array([occ[v][rng.integers(0,3)] for v in range(n)])
    Xh = np.bincount(alpha[om], minlength=g); mx = max(mx, (a - Xh).max())
print(f"n={n} s={s} g={g}: random omega: max_h(a-X) = {mx} = {mx/math.sqrt(b*s*ln(n)):.2f} sqrt(b s ln n); sqrt(6bs) = {math.sqrt(6*b*s):.1f}")
# adversarial omega for THIS point: avoid group 0
om = np.array([next((e for e in occ[v] if alpha[e] != 0), occ[v][0]) for v in range(n)])
Xh = np.bincount(alpha[om], minlength=g)
r = math.ceil(9*math.sqrt(b*s*ln(n)))
print(f"single-point adversarial omega: a_0={a[0]}, X_0={Xh[0]}, a_0-X_0={a[0]-Xh[0]} (r={r}; r>a_h at this n: {r>a.max()})")
# (5) Bernstein tail sanity at a moderate threshold: P(a_h - X_h >= lam) vs exp(-lam^2/(24 b s))
#     for a FIXED omega (worst-ish: fixed first occurrence), over fresh designs
om = np.array([occ[v][0] for v in range(n)])
dev = []
for _ in range(4000):
    lab = rng.integers(0, g, n); ue = cl[np.arange(n), rng.integers(0, 3, n)]; alpha = lab[ue]
    a = np.bincount(lab, minlength=g); Xh = np.bincount(alpha[om], minlength=g)
    dev.append(a[0]-Xh[0])
dev = np.array(dev)
print("fixed omega, var(a_0-X_0) =", round(dev.var(),1), " proof bound 6bs =", 6*b*s)
for c in (1, 2, 3):
    lam = c*math.sqrt(b*s)
    print(f"  P(a-X >= {c} sqrt(bs)) emp = {np.mean(dev>=lam):.4f}  Bernstein-as-used exp(-lam^2/(24bs)) = {math.exp(-lam*lam/(24*b*s)):.4f}")
