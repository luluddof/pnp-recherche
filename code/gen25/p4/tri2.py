# N=2: a family of 2x2 matrices is triangularisable by equivalence (over C) iff its span has dim <= 2, or dim 3 with a
# rank-1 annihilator Z (<F,Z>=0), i.e. iff some (r,c) != 0 has r^T F_m c = 0 for all m.
import sys, pickle, numpy as np
res = pickle.load(open(sys.argv[1], 'rb'))
for i, R in enumerate(res):
    out = []
    for side in ['Es', 'Cs']:
        F = R[side]; F = F/np.abs(F).max(); M = F.reshape(len(F), 4)
        u, s, vt = np.linalg.svd(M); s = s/s[0]
        k = (s > 1e-6).sum()
        info = "dim %d sv %s" % (k, np.round(s, 4))
        if k == 3:
            Z = vt[3].reshape(2, 2); zs = np.linalg.svd(Z, compute_uv=False)
            info += " annihilator sv ratio %.3f (%s)" % (zs[1]/zs[0], "NOT triangularisable" if zs[1]/zs[0] > 1e-4 else "triangularisable")
        elif k == 4: info += " NOT triangularisable"
        else: info += " triangularisable"
        out.append(side + ": " + info)
    print("sol", i, " | ".join(out))
