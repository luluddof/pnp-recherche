import time
from mysat import solve
from core import gstar
N = 5
for q in range(1, 11):
    p = 0
    while N*N > q*gstar(p, N) + p*gstar(q, N): p += 1
    lb = p
    while True:
        t = time.time(); r = solve(N, q, p, 0)
        if r: break
        p += 1
    print(f"N=5 ST q<={q}: min p = {p}  (counting bound {lb})", flush=True)
