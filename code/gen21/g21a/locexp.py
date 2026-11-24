# Local-exponent lemma (notes sect. 3): with S(E) = min_b [log Z(b) + b E],  d log S / d log E < 1/2  <=>  d/db [b log Z(b)] > 0.
# Scan of g(b) = d/db [b log Z(b)] per generator tower (modes r = D, D+1, ...), exact sums.
import math
def psi(kind, x):   # phi + x phi'
    if kind == 'f': return math.log1p(math.exp(-x)) - x/(math.exp(x)+1)
    return -math.log(-math.expm1(-x)) - x/math.expm1(x)
def g(kind, D, b, extra=0.0, K=20000):
    return extra + sum(psi(kind, b*(k+D)) for k in range(K) if b*(k+D) < 60)
bs = [0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 4, 6, 8, 12]
rows = [("fermion D=1/2 (NS)", 'f', 0.5, 0), ("fermion D=1 (sympl. fermion)", 'f', 1.0, 0), ("fermion Ramond: modes 1,2,.. + half zero mode", 'f', 1.0, 0.5*math.log(2)),
        ("boson D=1 (Heisenberg)", 'b', 1.0, 0), ("boson D=1/2 (sympl. boson / twisted Heis.)", 'b', 0.5, 0), ("fermion D=1/4", 'f', 0.25, 0), ("fermion D=3/4", 'f', 0.75, 0)]
print("g(b) = d/db[b log Z(b)] per tower;  >0 means local growth exponent of the entropy < 1/2 at that temperature")
print(" "*48 + " ".join(f"{b:>7}" for b in bs))
for name, kind, D, ex in rows:
    print(f"{name:48s}" + " ".join(f"{g(kind, D, b, ex):7.3f}" for b in bs))
print("pair fermion D=1/4 + D=3/4 (bc(1/4,3/4)):".ljust(48) + " ".join(f"{g('f',0.25,b)+g('f',0.75,b):7.3f}" for b in bs))
