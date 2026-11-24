# enumerate all counting-feasible (p,q) specs with max(p,q) <= 28 at N=16 (stabilisers of order 1,2,4, normal), run per group
import sys, itertools, multiprocessing as mp
from orbit import *; from groups import catalogue
N = 16
def pal_specs(maxc):
    out = []
    for n1 in range(0, 2):
        for n2 in range(0, 4):
            for n4 in range(0, 7):
                c = 16 * n1 + 8 * n2 + 4 * n4
                if 0 < c <= maxc: out.append((n1, n2, n4, c))
    return out
def specs(maxc=24):
    res = []
    for S in pal_specs(maxc):
        for T in pal_specs(maxc):
            if S > T: continue    # S <-> T symmetric
            cS = sum(g for g in [gstar(T[3])] * S[0]) + S[1] * (gstar(T[3]) // 2) + S[2] * (gstar(T[3]) // 4)
            cT = T[0] * gstar(S[3]) + T[1] * (gstar(S[3]) // 2) + T[2] * (gstar(S[3]) // 4)
            if cS + cT >= 16: res.append((S, T))
    return res
def child(g, sp_, caps, q): q.put(run(g, catalogue(N)[g], sp_, caps=caps, tries=300, verify=True))
if __name__ == "__main__":
    gl = sys.argv[1].split(','); to = int(sys.argv[2])
    SP = specs(); print("feasible palette specs (n_free, n_H2, n_H4, colours) S|T:", SP, flush=True)
    for g in gl:
        mul, inv = catalogue(N)[g][:2]; H1 = frozenset([0])
        H2 = normal_subgroups(mul, inv, 2)[:1]; H4 = normal_subgroups(mul, inv, 4)[:1]
        for S, T in SP:
            for h2 in H2 or [None]:
                for h4 in H4 or [None]:
                    if (S[1] or T[1]) and h2 is None: continue
                    if (S[2] or T[2]) and h4 is None: continue
                    sp_ = [('S', H1)] * S[0] + [('S', h2)] * S[1] + [('S', h4)] * S[2] + [('T', H1)] * T[0] + [('T', h2)] * T[1] + [('T', h4)] * T[2]
                    nT = T[3]; nS = S[3]
                    caps = [gstar(nT if p == 'S' else nS) // len(H) for p, H in sp_]
                    q = mp.Queue(); pr = mp.Process(target=child, args=(g, sp_, caps, q)); pr.start(); pr.join(to)
                    if pr.is_alive(): pr.kill(); pr.join(); print(f"  TIMEOUT {to}s ({g} S={S} T={T})", flush=True)
                    if not (S[1] or T[1]) and not (S[2] or T[2]): break
                if not (S[1] or T[1]): break
