import sys, time
from shape import *
N = int(sys.argv[1]); pmax = int(sys.argv[2]); tmo = float(sys.argv[3]) if len(sys.argv) > 3 else 1e9
# feasible(p,q) in the shape iff exists c: SAT(q-c, p-c, c).  By symmetry only p >= q.
t00 = time.time()
for q in range(1, pmax+1):
    for p in range(q, pmax+1):
        if time.time() - t00 > tmo: print("time out"); sys.exit()
        res = None
        for c in range(0, q+1):
            t0 = time.time(); cells = sat_shape(N, q-c, p-c, c)
            if cells:
                ok, pp, qq = verify_shape(N, cells, quiet=True); assert ok and pp <= p and qq <= q
                res = (c, pp, qq); break
        print(f"N={N} (p,q)=({p},{q}): " + (f"FEASIBLE in shape with U={res[0]} (exact ranks {res[1]},{res[2]})" if res else "infeasible in shape (all U)"), f"[{time.time()-t00:.0f}s]", flush=True)
        if res: break
