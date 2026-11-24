import sys, time
from interp import test
from monos import monos
W = [(5,(3,1,1),(2,2,1)),(5,(2,2,1),(2,2,1)),(6,(4,1,1),(4,1,1)),(6,(3,2,1),(4,1,1)),(6,(4,1,1),(2,2,2)),(6,(3,2,1),(3,2,1))]
for d, mu, nu in W:
    cont = (mu,nu,nu,mu); t0 = time.time()
    m, a, b = test(cont)
    print(d, cont, f"monomials={m} eqs_of_J={a} nonzero_on_Z3={b}", 'WITNESS' if b > 0 else '', f"{time.time()-t0:.0f}s", flush=True)
