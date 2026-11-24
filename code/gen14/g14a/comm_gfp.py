"""G14a: commutative class over GF(p) (p odd), exact min total #bad edges T_p(l,d;tau) for EVERY tau-set
(up to scaling t->lambda t and permutations), by lazy hitting-set MILP (as G8d/lazy.py) with feasibility
checked exactly mod p.  P_z(t)=1+sum u_i(z)t^i, edge (z,z^k) good iff (P_z+P_{z^k})(tau_k)=0.
usage: python3 comm_gfp.py l d p [timelimit_per_class]"""
import sys, itertools, time
import numpy as np, flint
from scipy.optimize import milp, LinearConstraint, Bounds

def edges(l):
    return [(z, z | (1 << k), k) for z in range(1 << l) for k in range(l) if not z >> k & 1]

class Inst:
    def __init__(s, l, d, p, tau):
        s.l, s.d, s.p, s.tau = l, d, p, tau
        s.E = edges(l); s.idx = {(z, k): i for i, (z, w, k) in enumerate(s.E)}
        s.n = (1 << l) * d
    def rows(s, es):
        d = s.d; P = s.p; out = []
        for e in es:
            z, w, k = s.E[e]; r = [0] * (s.n + 1); tk = s.tau[k]; pw = tk
            for j in range(d):
                r[z * d + j] = pw; r[w * d + j] = pw; pw = pw * tk % P
            r[s.n] = (P - 2) % P; out.append(r)
        return out
    def feasible(s, es):
        if not es: return True
        R = s.rows(es)
        A = flint.nmod_mat([r[:s.n] for r in R], s.p); Ab = flint.nmod_mat(R, s.p)
        return A.rank() == Ab.rank()
    def certificate(s, es):
        cur = list(es); ch = max(1, len(cur) // 4)
        while ch >= 1:
            i = 0
            while i < len(cur):
                trial = cur[:i] + cur[i + ch:]
                if not s.feasible(trial): cur = trial
                else: i += ch
            ch //= 2
        return cur
    def images(s, cert):   # translations only (valid for every tau)
        out = set()
        for c in range(1 << s.l):
            out.add(tuple(sorted(s.idx[((s.E[e][0] ^ c) & ~(1 << s.E[e][2]), s.E[e][2])] for e in cert)))
        return out

def solve(l, d, p, tau, tl=300):
    I = Inst(l, d, p, tau); ne = len(I.E); cons = set(); t0 = time.time(); it = 0
    for S in itertools.combinations(range(l), min(d + 1, l)):
        cert = []
        for sub in range(1 << len(S)):
            z = sum(1 << S[i] for i in range(len(S)) if sub >> i & 1)
            for k in S:
                if not z >> k & 1: cert.append(I.idx[(z, k)])
        if not I.feasible(cert): cons |= I.images(cert)
    while True:
        it += 1
        if time.time() - t0 > tl: return None, lb if it > 1 else 0, it
        C = sorted(cons)
        if C:
            A = np.zeros((len(C), ne))
            for i, c in enumerate(C): A[i, list(c)] = 1
            res = milp(np.ones(ne), constraints=LinearConstraint(A, 1, np.inf), integrality=np.ones(ne),
                       bounds=Bounds(0, 1), options={"time_limit": 120})
            if res.x is None: return None, 0, it
            x = res.x
        else: x = np.zeros(ne)
        bad = [e for e in range(ne) if x[e] > .5]; good = [e for e in range(ne) if x[e] <= .5]
        lb = len(bad)
        if I.feasible(good):
            per = [sum(1 for e in bad if I.E[e][2] == k) for k in range(l)]
            return len(bad), per, it
        cons |= I.images(I.certificate(good))

def classes(l, p):
    seen = set(); out = []
    for S in itertools.combinations(range(1, p), l):
        if S in seen: continue
        orb = {tuple(sorted(lam * s % p for s in S)) for lam in range(1, p)}
        seen |= orb; out.append(S)
    return out

if __name__ == "__main__":
    l, d, p = map(int, sys.argv[1:4]); tl = float(sys.argv[4]) if len(sys.argv) > 4 else 300
    best = None; cl = classes(l, p)
    print(f"(l,d,p)=({l},{d},{p}): {len(cl)} tau classes", flush=True)
    for tau in cl:
        t0 = time.time(); T, per, it = solve(l, d, p, list(tau), tl)
        if T is None:
            print(f"  tau={tau}: TIMEOUT, relaxation LB {per}, iters {it}", flush=True); continue
        print(f"  tau={tau}: T={T} per={per} iters={it} {time.time()-t0:.1f}s", flush=True)
        if best is None or T < best[0]: best = (T, tau, per)
    print(f"MIN over tau classes (l,d,p)=({l},{d},{p}): {best}", flush=True)
