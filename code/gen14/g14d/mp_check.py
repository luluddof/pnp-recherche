# G14d: multi-point bounds, exact small-case optima in the commutative class, adversary optimisation.
import itertools, sys
from fractions import Fraction as Fr
from math import comb
P = (1<<61)-1

def bound_old(l, N, d, mc):   # G11b: sum_x max(0, m+1-mu(x)), m=floor((N-d-1)/2)
    if N-d-1 < 0: return 0
    m = (N-d-1)//2
    return sum(max(0, m+1-sum(mc[c] for c in range(l) if x>>c&1)) for x in range(1<<l))
def bound_new(l, N, d, mc):   # G14d: sum_x max(0, N-d-2mu(x))
    return sum(max(0, N-d-2*sum(mc[c] for c in range(l) if x>>c&1)) for x in range(1<<l))

def rank_consistent(rows, rhs, ncol, mod=P):
    """return True iff system rows*x = rhs consistent (mod p)."""
    M = [r[:] + [b] for r, b in zip(rows, rhs)]
    piv = 0
    for col in range(ncol):
        pr = None
        for i in range(piv, len(M)):
            if M[i][col] % mod: pr = i; break
        if pr is None: continue
        M[piv], M[pr] = M[pr], M[piv]
        inv = pow(M[piv][col], mod-2, mod)
        M[piv] = [(v*inv) % mod for v in M[piv]]
        for i in range(len(M)):
            if i != piv and M[i][col] % mod:
                f = M[i][col]
                M[i] = [(a - f*b) % mod for a, b in zip(M[i], M[piv])]
        piv += 1
        if piv == len(M): break
    for i in range(piv, len(M)):
        if M[i][ncol] % mod: return False
    return True

def consistent_exact(rows, rhs, ncol):
    M = [[Fr(v) for v in r] + [Fr(b)] for r, b in zip(rows, rhs)]
    piv = 0
    for col in range(ncol):
        pr = None
        for i in range(piv, len(M)):
            if M[i][col] != 0: pr = i; break
        if pr is None: continue
        M[piv], M[pr] = M[pr], M[piv]
        inv = 1/M[piv][col]
        M[piv] = [v*inv for v in M[piv]]
        for i in range(len(M)):
            if i != piv and M[i][col] != 0:
                f = M[i][col]
                M[i] = [a - f*b for a, b in zip(M[i], M[piv])]
        piv += 1
        if piv == len(M): break
    return all(M[i][ncol] == 0 for i in range(piv, len(M)))

def chi(z): return -1 if bin(z).count('1') & 1 else 1

def equations(l, d, pts):
    """pairs (edge z (z_c=0), point k); equation: sum_j (G[z][j]-G[z^c][j]) tau^j = -(chi(z)-chi(z^c)) = -2chi(z)"""
    eqs = []
    for k, (tau, c) in enumerate(pts):
        for z in range(1<<l):
            if z>>c&1: continue
            row = [0]*((1<<l)*d)
            for j in range(1, d+1):
                row[z*d+j-1] += tau**j
                row[(z^(1<<c))*d+j-1] -= tau**j
            eqs.append((row, -2*chi(z), k))
    return eqs

def comm_opt(l, d, pts, maxb=None, exact_all=False):
    eqs = equations(l, d, pts); n = len(eqs); ncol = (1<<l)*d
    if maxb is None: maxb = n
    chk = consistent_exact if exact_all else rank_consistent
    for b in range(0, maxb+1):
        for bad in itertools.combinations(range(n), b):
            bs = set(bad)
            rows = [eqs[i][0] for i in range(n) if i not in bs]
            rhs = [eqs[i][1] for i in range(n) if i not in bs]
            if d == 0:
                ok = (len(rows) == 0)
            else:
                ok = chk(rows, rhs, ncol)
            if ok:
                if d > 0: assert consistent_exact(rows, rhs, ncol)   # exact witness
                prof = [sum(1 for i in bad if eqs[i][2] == k) for k in range(len(pts))]
                return b, prof
    return None, None

if __name__ == "__main__":
    print("== A. exact commutative optimum T vs bounds (old=G11b, new=G14d). 'T' exact (infeasibility of smaller bad sets mod 2^61-1, l=1,2 exact over Q)")
    cases = []
    for l in (1, 2):
        for N in range(1, 7 if l == 1 else 6):
            for mc in itertools.product(range(N+1), repeat=l):
                if sum(mc) != N or list(mc) != sorted(mc, reverse=True): continue
                for d in range(0, N):
                    cases.append((l, N, d, mc))
    for l, N, d, mc in [(3,3,1,(1,1,1)), (3,3,2,(1,1,1)), (3,4,1,(2,1,1)), (3,4,2,(2,1,1)), (3,4,2,(2,2,0)), (3,4,3,(2,1,1)),
                        (3,4,1,(4,0,0)), (3,5,3,(2,2,1)), (3,5,4,(2,2,1)), (3,5,2,(3,1,1)), (3,4,1,(2,2,0)), (3,5,3,(5,0,0)), (3,5,3,(3,2,0)), (3,6,4,(2,2,2)), (3,6,5,(2,2,2))]:
        cases.append((l, N, d, mc))
    allok = True
    for l, N, d, mc in cases:
        pts = []; t = 1
        for c in range(l):
            for _ in range(mc[c]): pts.append((t, c)); t += 1
        bo, bn = bound_old(l, N, d, mc), bound_new(l, N, d, mc)
        T, prof = comm_opt(l, d, pts, exact_all=(l <= 2), maxb=11)
        ok = (T is None) or (T >= bn and T >= bo)
        allok &= ok
        print(f"l={l} N={N} d={d} m_c={mc}: old={bo} new={bn} T_comm={T if T is not None else '>11'} profile={prof} {'OK' if ok else 'VIOLATION'}")
        sys.stdout.flush()
    print("ALL OK" if allok else "SOME VIOLATION")
