#!/usr/bin/env python3
"""G18S r11: 'u=2 full packet' GF(2) search, profiles (3,2) at N=5 and (3,3),(4,2),(5,1) at N=6, 4 in parallel, cap argv[1] s."""
import sys, time, multiprocessing as mp
from dflib import *
CAP = int(sys.argv[1]); u = 2; wires = [(0, 3), (3, 0), (1, 2), (2, 1)]
def worker(d, q):
    t = time.time(); r, tabs = sat_search(u, wires, d)
    q.put((d, r, verify_gf2(u, wires, d, tabs) if r == 'SAT' else None, round(time.time() - t, 1)))
if __name__ == "__main__":
    q = mp.Queue(); ps = []
    for d in [(0, 0, 0, 1, 1), (0, 0, 0, 1, 1, 1), (0, 0, 0, 0, 1, 1), (0, 0, 0, 0, 0, 1)]:
        p = mp.Process(target=worker, args=(d, q)); p.start(); ps.append((p, d))
    t0 = time.time(); got = set()
    while time.time() - t0 < CAP and len(got) < 4:
        time.sleep(1)
        while not q.empty():
            d, r, ok, tt = q.get(); got.add(d); print(f"N={len(d)} dirs={d}: {r} {'(verified)' if ok else ''} {tt}s", flush=True)
    for p, d in ps:
        if p.is_alive(): p.terminate(); print(f"N={len(d)} dirs={d}: TIMEOUT {CAP}s", flush=True)
