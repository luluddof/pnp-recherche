# driver: python drv.py N mode timeout group1,group2,...   (each spec in a child process with a wall-clock timeout)
import sys, itertools, multiprocessing as mp, time
from orbit import *; from groups import catalogue
def child(N, gname, sp_, caps, q):
    G = catalogue(N)[gname]
    q.put(run(gname, G, sp_, caps=caps, tries=300))
def specs_for(N, mul, inv, mode):
    H1 = frozenset([0]); H2s = normal_subgroups(mul, inv, 2); H4s = normal_subgroups(mul, inv, 4)
    pairs = lambda L: [(a, b) for i, a in enumerate(L) for b in L[i:]][:12]
    out = []
    if mode.startswith('free'):       # free_nS_nT
        _, nS, nT = mode.split('_'); out.append([('S', H1)] * int(nS) + [('T', H1)] * int(nT))
    elif mode == 'f+h2':
        for HS, HT in pairs(H2s): out.append([('S', H1), ('S', HS), ('T', H1), ('T', HT)])
    elif mode == 'f+h4h4':
        for HS, HT in pairs(H4s): out.append([('S', H1), ('S', HS), ('S', HS), ('T', H1), ('T', HT), ('T', HT)])
    elif mode == 'h2x3':
        for HS, HT in pairs(H2s): out.append([('S', HS)] * 3 + [('T', HT)] * 3)
    elif mode == 'h2x3f':   # S: free + H2 ; T: H2 x 3  -> (24,24) at 16
        for HS in H2s:
            for HT in H2s: out.append([('S', H1), ('S', HS), ('T', HT), ('T', HT), ('T', HT)])
    return out
if __name__ == "__main__":
    N = int(sys.argv[1]); mode = sys.argv[2]; to = int(sys.argv[3]); gl = sys.argv[4].split(',')
    for gname in gl:
        G = catalogue(N)[gname]; mul, inv = G[0], G[1]
        for sp_ in specs_for(N, mul, inv, mode):
            nT = sum(N // len(H) for p, H in sp_ if p == 'T'); nS = sum(N // len(H) for p, H in sp_ if p == 'S')
            caps = [gstar(nT if p == 'S' else nS) // len(H) for p, H in sp_]
            q = mp.Queue(); pr = mp.Process(target=child, args=(N, gname, sp_, caps, q)); pr.start(); pr.join(to)
            if pr.is_alive():
                pr.kill(); pr.join(); print(f"  TIMEOUT {to}s  ({gname} {[(p, sorted(H)) for p, H in sp_]})", flush=True)
