# SU(3,q) Sylow ansatz with t(a, *) in {(mu a, *)}: mod-Z part linear.  python ztmu.py q timeout
import sys, multiprocessing as mp
from zt import *; from g64 import su3_sylow, gf
q = int(sys.argv[1]); to = int(sys.argv[2])
G = su3_sylow(q); mul, inv, el = G; N = len(mul)
n = {2: 2, 4: 4, 8: 6}[q]; poly = {2: 0b111, 4: 0b10011, 8: 0b1000011}[q]; m = gf(n, poly)
Z = frozenset(x for x in range(N) if el[x][0] == 0)
Zs = sorted(Z); Z0 = None
for cand in central_subgroups(mul, inv, len(Z) // 2):
    if cand <= Z: Z0 = cand; break
def child(mu):
    E0 = ZT.__new__(ZT)
    # cosets as in ZT: compute allowed per coset index after building cosets the same way
    cos = {}; reps = []
    for x in range(N):
        k = min(mul[x][z] for z in Zs)
        if k not in cos: cos[k] = len(reps); reps.append(k)
    allowed = [[x for x in range(N) if el[x][0] == m(mu, el[reps[c]][0])] for c in range(len(reps))]
    runzt(f"SU3({q})_mu{mu}", G, Z, Z0, allowed=allowed)
mus = [int(a) for a in sys.argv[3].split(',')] if len(sys.argv) > 3 else range(2, 1 << n)
for mu in mus:
    pr = mp.Process(target=child, args=(mu,)); pr.start(); pr.join(to)
    if pr.is_alive(): pr.kill(); pr.join(); print(f"  TIMEOUT {to}s mu={mu}", flush=True)
