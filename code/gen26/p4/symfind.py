# support-level automorphisms of a MC solution: (A,B) perms, term perm sigma, flip f: E_l -> A E_l B^T = E_sig(l) etc.
import sys, json, itertools, numpy as np
d = json.load(open(sys.argv[1])); kE, kC = sys.argv[2], sys.argv[3]
Es = [np.abs(np.array(m, dtype=float).astype(int)) for m in d[kE]]; Cs = [np.abs(np.array(m, dtype=float).astype(int)) for m in d[kC]]; N = Es[0].shape[0]; q = len(Es)
keyE = {Es[l].tobytes() + Cs[l].tobytes(): l for l in range(q)}
P = [np.eye(N, dtype=Es[0].dtype)[list(pm)] for pm in itertools.permutations(range(N))]
auts = []
for flip in (0, 1):
    for A in P:
        AE = [A @ (Cs[l] if flip else Es[l]) for l in range(q)]
        for B in P:
            sig = []
            for l in range(q):
                k = (AE[l] @ B.T).tobytes() + (B @ (Es[l] if flip else Cs[l]) @ A.T).tobytes()
                if k not in keyE: break
                sig.append(keyE[k])
            else:
                auts.append((flip, tuple(np.argmax(A, 1)), tuple(np.argmax(B, 1)), tuple(sig)))
print(len(auts), "support automorphisms (incl identity)")
for a in auts[:20]: print(a)
