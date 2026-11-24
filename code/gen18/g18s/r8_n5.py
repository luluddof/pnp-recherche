#!/usr/bin/env python3
"""G18S r8: GF(2) search at N = 5, 6 for 'u=2 full packet' and 'u=3 four pairs', one SAT call per direction profile,
4 in parallel, wall-clock cap per call (argv[1])."""
import itertools, sys, time, multiprocessing as mp
from dflib import *
CAP = int(sys.argv[1])
def job(args):
    name, u, wires, dirs = args
    t = time.time(); r, tabs = sat_search(u, wires, dirs)
    return name, dirs, r, (verify_gf2(u, wires, dirs, tabs) if r == 'SAT' else None), round(time.time() - t, 1)
if __name__ == "__main__":
    full = ("u=2 full packet", 2, [(0, 3), (3, 0), (1, 2), (2, 1)]); four = ("u=3 four pairs", 3, [(0, 7), (1, 6), (2, 5), (3, 4)])
    jobs = []
    for N in (5, 6):
        jobs += [(full[0], full[1], full[2], d) for d in itertools.combinations_with_replacement(range(2), N) if d.count(0) >= d.count(1)]
    jobs += [(four[0], four[1], four[2], d) for d in [(0, 0, 0, 1, 2), (0, 1, 1, 1, 2), (0, 1, 2, 2, 2), (0, 0, 1, 1, 2), (0, 0, 0, 1, 1), (0, 0, 0, 0, 1)]]
    procs = []; pending = list(jobs); running = []
    q = mp.Queue()
    def worker(a, q): q.put(job(a))
    t0 = time.time()
    while pending or running:
        while pending and len(running) < 4:
            a = pending.pop(0); p = mp.Process(target=worker, args=(a, q)); p.start(); running.append((p, a, time.time()))
        time.sleep(0.5)
        while not q.empty():
            name, dirs, r, ok, tt = q.get(); print(f"{name} N={len(dirs)} dirs={dirs}: {r} {'(verified)' if ok else ''} {tt}s", flush=True)
        for (p, a, ts) in list(running):
            if not p.is_alive(): running.remove((p, a, ts))
            elif time.time() - ts > CAP:
                p.terminate(); running.remove((p, a, ts)); print(f"{a[0]} N={len(a[3])} dirs={a[3]}: TIMEOUT {CAP}s", flush=True)
    while not q.empty():
        name, dirs, r, ok, tt = q.get(); print(f"{name} N={len(dirs)} dirs={dirs}: {r} {'(verified)' if ok else ''} {tt}s", flush=True)
