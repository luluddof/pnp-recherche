import sys, time, json
sys.path.insert(0,'/home/user/pnp-recherche/code/gen22/g22b')
from shape import sat_shape, verify_shape
N,q,p=map(int,sys.argv[1:4]); ms=int(sys.argv[4]) if len(sys.argv)>4 else None
t=time.time(); cells=sat_shape(N,q,p,0,maxsize=ms)
print(f"ST N={N} S(q)={q} T(p)={p} maxsize={ms}:", "SAT" if cells else "UNSAT", f"{time.time()-t:.0f}s",flush=True)
if cells:
    ok,pp,qq=verify_shape(N,cells); print("verify",ok,pp,qq)
    json.dump({"N":N,"p":pp,"q":qq,"cells":cells},open(f"st_{N}_{pp}_{qq}.json","w"))
