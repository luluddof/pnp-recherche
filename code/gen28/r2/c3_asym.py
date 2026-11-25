import numpy as np
def kl(a, b):
    a = np.clip(a, 1e-300, 1-1e-16)
    return a*np.log(a/b) + (1-a)*np.log((1-a)/(1-b))
for n in [2, 3]:
    b = (n-1)/(2*n); c = np.log(n)/np.log(n*n-1)
    C = n*np.log(n)/(2*(n-1)*np.log(n*n-1)**2)
    for eta in [0.05, 0.1, 0.2, 0.3, 0.5]:
        ys = np.linspace(1e-4, 1, 200001); Dk = c*eta
        ok = ys > Dk
        rate = np.where(ok, kl(ys, 2*b) + ys*kl(0.5 + Dk/(2*np.maximum(ys,1e-9)), 0.5), np.inf)
        dl = rate.min()/np.log(n)
        print(f"n={n} eta={eta}: exact-limit delta bound {dl:.5f}  vs C*eta^2 = {C*eta**2:.5f} (C={C:.3f}), ratio {dl/(C*eta**2):.3f}")
