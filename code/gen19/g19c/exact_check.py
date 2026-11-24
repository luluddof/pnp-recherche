"""Exact (integer) verification of the product point for n=(n01,n02,n12): P = -Y(x)Y(x)Y + star0 + star1 + star2, Y = 1 - SWAP."""
import numpy as np, sys
from starlib import *
n = tuple(int(a) for a in sys.argv[1].split(',')); inst = Inst(*n)
def swap(m):
    S = np.zeros((m, m, m, m), dtype=np.int64)
    for a in range(m):
        for b in range(m): S[b, a, a, b] = 1
    return S
def ident(m): return np.einsum('ac,bd->abcd', np.eye(m, dtype=np.int64), np.eye(m, dtype=np.int64))
def asm(F01, F02, F12): return np.einsum('uvab,UwAc,VWBC->uUvVwWaAbBcC', F01, F02, F12).reshape(inst.D, inst.D)
I = [ident(m) for m in n]; S = [swap(m) for m in n]; Y = [i - s for i, s in zip(I, S)]
Pint = asm(*S); print('P = S(x)S(x)S equals inst.P:', np.array_equal(Pint, inst.P.astype(np.int64)))
A = -asm(Y[0], Y[1], Y[2])
st0 = asm(Y[0], Y[1], I[2]) - asm(Y[0], I[1], I[2]) + asm(I[0], I[1], I[2])
st1 = asm(Y[0], I[1], Y[2]) - asm(I[0], I[1], Y[2])
st2 = asm(I[0], Y[1], Y[2]) - asm(I[0], Y[1], I[2])
print('exact integer identity P = A + st0 + st1 + st2:', np.array_equal(Pint, A + st0 + st1 + st2))
print('rank A =', np.linalg.matrix_rank(A.astype(float)), ' product value', inst.prod_r, ' |A|_F^2 =', int((A*A).sum()))
xs = []
for s, st in enumerate([st0, st1, st2]):
    c, l, m, ba, bb = inst.stars[s]; T3 = inst.to3(st.astype(float))
    mr = [np.linalg.matrix_rank(np.moveaxis(T3, k, 0).reshape(inst.K[k], -1)) for k in range(3)]
    print(f' star {s} (centre {c}): multilinear ranks {mr}; allowed bonds at leaves {l},{m}: {ba},{bb}')
    G = np.ascontiguousarray(inst._unarr(c, T3))   # (c,l,m)
    LA = np.linalg.svd(G.transpose(1, 0, 2).reshape(inst.K[l], -1), full_matrices=False)[0][:, :ba]
    LB = np.linalg.svd(G.transpose(2, 0, 1).reshape(inst.K[m], -1), full_matrices=False)[0][:, :bb]
    C = np.einsum('clm,la,mb->cab', G, LA, LB); xs += [C.ravel(), LA.ravel(), LB.ravel()]
x = np.concatenate(xs); nm = norms(inst, x, inst.prod_r)
print(' star-parameter form reproduces: tail =', nm['tail'], 'norms', nm)
np.savez(f"product_n{''.join(map(str,n))}.npz", x=x)
