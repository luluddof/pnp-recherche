# python adrv.py N nS nT autoorder timeout groups [zpair]
import sys, multiprocessing as mp
from orbit import *; from groups import catalogue; from auto import automorphisms_of_order
def child(N, g, blocks, caps, phis, z):
    run(g + f'_aut{len(phis)}', catalogue(N)[g], blocks, caps=caps, tries=300, autos=phis, zpair=z)
if __name__ == "__main__":
    N, nS, nT, k, to = map(int, sys.argv[1:6]); gl = sys.argv[6].split(','); zp = len(sys.argv) > 7
    H1 = frozenset([0])
    blocks = ([('S', H1), ('T', H1)] * nS) if zp else ([('S', H1)] * nS + [('T', H1)] * nT)
    caps = [gstar(N * (nT if p == 'S' else nS)) for p, H in blocks]
    for g in gl:
        mul, inv = catalogue(N)[g][:2]
        for seed in range(3):
            A = automorphisms_of_order(mul, k, limit=1, seed=seed)
            if not A: print(g, 'no automorphism of order', k); break
            z = None
            if zp:
                z = next(h for h in range(1, N) if mul[h][h] == 0 and all(mul[h][x] == mul[x][h] for x in range(N)))
            print(g, 'automorphism', A[0], flush=True)
            pr = mp.Process(target=child, args=(N, g, blocks, caps, A, z)); pr.start(); pr.join(to)
            if pr.is_alive(): pr.kill(); pr.join(); print(f"  TIMEOUT {to}s ({g} seed {seed})", flush=True)
