import math, numpy as np
from opt import taylor, best_V3
print("V3 only, large n.  log2n  E/n/(L2/sqrtL) any-m (r,Q,m) | m=3^a")
for lg in [1e12,1e15,1e20,1e30,1e50,1e100,1e200]:
    L=lg*math.log(2); L2=math.log(L); u=L2/math.sqrt(L)
    b=best_V3(L); b3=best_V3(L,True)
    print(f"{lg:8.0e} {b[0]/u:8.3f} {b[1:]} | {b3[0]/u:8.3f} m={b3[3]}  m/sqrt(Q/lnQ)={b[3]/math.sqrt(b[2]/math.log(b[2])):.2f}",flush=True)
