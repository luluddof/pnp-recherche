# python zdrv.py N npairs timeout groups [exact_size]
import sys, multiprocessing as mp
from orbit import *; from groups import catalogue
def child(N, g, z, blocks, caps, exact): run(g + f'_z{z}', catalogue(N)[g], blocks, caps=caps, exact=exact, tries=300, zpair=z)
if __name__ == "__main__":
    N, npairs, to = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); gl = sys.argv[4].split(',')
    ex = int(sys.argv[5]) if len(sys.argv) > 5 else None
    H1 = frozenset([0]); blocks = [('S', H1), ('T', H1)] * npairs
    caps = [gstar(N * npairs)] * (2 * npairs); exact = [ex] * (2 * npairs) if ex else None
    for g in gl:
        mul, inv = catalogue(N)[g][:2]
        zs = [h for h in range(1, N) if mul[h][h] == 0 and all(mul[h][x] == mul[x][h] for x in range(N))]
        for z in zs[:3]:
            pr = mp.Process(target=child, args=(N, g, z, blocks, caps, exact)); pr.start(); pr.join(to)
            if pr.is_alive(): pr.kill(); pr.join(); print(f"  TIMEOUT {to}s ({g} z={z})", flush=True)
