"""G8d: exact minimum number of bad edges T(l,d) (generic tau) via lazy certificate generation.
Problem (commutative class, one point per direction): P_z(t)=1+sum_{i=1..d} u_i(z) t^i;
edge (z,z^k) good iff (P_z+P_{z^k})(tau_k)=0.  Good set G feasible iff no y on G with y^T A=0, y^T b!=0
(Fredholm).  Every infeasible support (certificate) must contain a bad edge.  We solve the hitting-set
MILP over the certificates found so far, test the optimum for algebraic feasibility (mod p, random tau =
generic tau w.h.p.), and if infeasible extract a minimal certificate (deletion filter), add all its images
under Aut(Q_l) (valid for generic tau), and repeat.  Final answer = optimum of a relaxation that is feasible
=> exact T (generic tau).  mode 'total' or 'max' (minimise max_k bad_k).
usage: python3 lazy.py l d [total|max]
"""
import sys, itertools, random, time
import numpy as np, flint
from scipy.optimize import milp, LinearConstraint, Bounds
P = 2**31 - 1

def edges(l):
    return [(z, z | (1 << k), k) for z in range(1 << l) for k in range(l) if not z >> k & 1]

class Inst:
    def __init__(s, l, d, seed=7):
        s.l, s.d = l, d
        rng = random.Random(seed)
        s.tau = [rng.randrange(2, P - 1) for _ in range(l)]
        s.E = edges(l); s.idx = {(z, k): i for i, (z, w, k) in enumerate(s.E)}
        s.n = (1 << l) * d
    def rows(s, es):
        d = s.d; out = []
        for e in es:
            z, w, k = s.E[e]; r = [0] * (s.n + 1); tk = s.tau[k]; pw = tk
            for j in range(d):
                r[z * d + j] = pw; r[w * d + j] = pw; pw = pw * tk % P
            r[s.n] = P - 2; out.append(r)
        return out
    def feasible(s, es):
        if not es: return True
        R = s.rows(es)
        A = flint.nmod_mat([r[:s.n] for r in R], P); Ab = flint.nmod_mat(R, P)
        return A.rank() == Ab.rank()
    def certificate(s, es):
        """minimal infeasible subset of es (assumed infeasible)."""
        cur = list(es)
        # shrink quickly: halving chunks, then single deletions
        ch = max(1, len(cur) // 4)
        while ch >= 1:
            i = 0
            while i < len(cur):
                trial = cur[:i] + cur[i + ch:]
                if not s.feasible(trial): cur = trial
                else: i += ch
            ch //= 2
        return cur
    def images(s, cert):
        l = s.l; out = set()
        for perm in itertools.permutations(range(l)):
            for c in range(1 << l):
                img = []
                for e in cert:
                    z, w, k = s.E[e]
                    zz = sum(1 << perm[i] for i in range(l) if z >> i & 1) ^ c
                    kk = perm[k]; base = zz & ~(1 << kk)
                    img.append(s.idx[(base, kk)])
                out.add(tuple(sorted(img)))
        return out

def solve(l, d, mode="total", tl=600, log=print):
    I = Inst(l, d); ne = len(I.E)
    cons = set()
    # seed: all (d+1)-subcubes (Barrier B certificates)
    for S in itertools.combinations(range(l), d + 1):
        cert = []
        for sub in range(1 << (d + 1)):
            z = sum(1 << S[i] for i in range(d + 1) if sub >> i & 1)
            for k in S:
                if not z >> k & 1: cert.append(I.idx[(z, k)])
        cons |= I.images(cert); break
    it = 0; t0 = time.time(); certs = []
    while True:
        it += 1
        C = sorted(cons); nv = ne + (1 if mode == "max" else 0)
        A = np.zeros((len(C) + (l if mode == "max" else 0), nv))
        for i, c in enumerate(C): A[i, list(c)] = 1
        lo = [1] * len(C); hi = [np.inf] * len(C); cost = np.zeros(nv)
        if mode == "max":
            for k in range(l):
                r = len(C) + k
                for e, (z, w, kk) in enumerate(I.E):
                    if kk == k: A[r, e] = 1
                A[r, -1] = -1; lo.append(-np.inf); hi.append(0)
            cost[-1] = 1; cost[:ne] = 1e-3   # tie-break: few total
        else: cost[:] = 1
        ub = np.ones(nv)
        if mode == "max": ub[-1] = np.inf
        res = milp(cost, constraints=LinearConstraint(A, lo, hi), integrality=np.ones(nv),
                   bounds=Bounds(0, ub), options={"time_limit": tl})
        bad = [e for e in range(ne) if res.x[e] > .5]
        good = [e for e in range(ne) if res.x[e] <= .5]
        per = [sum(1 for e in bad if I.E[e][2] == k) for k in range(l)]
        opt = res.status == 0
        if I.feasible(good):
            log(f"(l,d)=({l},{d}) {mode}: T={len(bad)} max_k={max(per)} per={per} iters={it} "
                f"#certs={len(certs)} rows={len(C)} milp_optimal={opt} time={time.time()-t0:.0f}s")
            return len(bad), per, bad, certs, opt
        cert = I.certificate(good); certs.append(cert)
        cons |= I.images(cert)
        if it % 5 == 0:
            log(f"  it {it}: relax bound {len(bad)} (max {max(per)}), cert size {len(cert)}, rows {len(cons)}, "
                f"{time.time()-t0:.0f}s", flush=True) if False else print(
                f"  it {it}: relax {len(bad)} (max {max(per)}), cert |supp|={len(cert)}, rows {len(cons)}, "
                f"{time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    l, d = int(sys.argv[1]), int(sys.argv[2]); mode = sys.argv[3] if len(sys.argv) > 3 else "total"
    T, per, bad, certs, opt = solve(l, d, mode)
    from collections import Counter
    print("certificate sizes:", sorted(Counter(len(c) for c in certs).items()))
    I = Inst(l, d)
    for c in certs[:6]:
        print("  cert:", sorted((bin(I.E[e][0])[2:].zfill(l), I.E[e][2]) for e in c))
