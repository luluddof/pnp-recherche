# check the anisotropy criterion (squares of non-central elements avoid Z0) against ZT SAT outcomes, order 8 and 16
from groups import catalogue
from zt import central_subgroups, ZT
from pysat.solvers import Cadical195
for N in (8, 16):
    for g, G in catalogue(N).items():
        mul, inv = G[0], G[1]; m = 2 if N == 8 else 4
        for Z in central_subgroups(mul, inv, m)[:4]:
            quot_elab = all(mul[x][x] in Z for x in range(N))
            for Z0 in central_subgroups(mul, inv, m // 2):
                if not Z0 <= Z: continue
                anis = all(mul[x][x] not in Z0 for x in range(N) if x not in Z)
                E = ZT(mul, inv, Z, Z0); S = Cadical195(bootstrap_with=E.cl); r = S.solve(); S.delete()
                if r or anis: print(N, g, sorted(Z), sorted(Z0), 'Gamma/Z elem.ab:', quot_elab, 'anisotropic:', anis, 'ZT SAT:', r)
print('(rows printed only where SAT or anisotropic; all other (Z,Z0) are non-anisotropic and UNSAT)')
