"""G18d: pruning lemma check. Any feasible B contains a feasible row-basis B' (|B'| = rank ev_B on the dual space);
on each line L_z at most l-d nodes are independent.  Build B' from the layered set, round-robin over directions
(balanced), verify by the PRIMAL system mod p (independent of the dual formulation)."""
import sys, itertools
from quad import *
def incr_basis(l, d, tau, order):
    Dp = l - d - 1; bs = basis(l, Dp); inv = [pow(t % P, P - 2, P) for t in tau]
    piv = {}; chosen = []
    for z, k in order:
        r = row_mod(l, tau, inv, bs, z, k)
        for c in sorted(piv):
            if r[c]:
                f = r[c]; pr = piv[c]; r = [(a - f * b) % P for a, b in zip(r, pr)]
        nz = next((i for i, a in enumerate(r) if a), None)
        if nz is None: continue
        iv = pow(r[nz], P - 2, P); r = [a * iv % P for a in r]
        piv[nz] = r; chosen.append((z, k))
    return chosen
for a in sys.argv[1:]:
    l, d = map(int, a.split(",")); tau = taus_for(l, "rand")
    L = layered(l, d)
    byk = [[e for e in L if e[1] == k] for k in range(l)]
    order = [e for grp in itertools.zip_longest(*byk) for e in grp if e is not None]
    Bp = incr_basis(l, d, tau, order)
    rk, ok = exact_dual(l, d, tau, Bp); pr = exact_primal(l, d, tau, Bp) if l <= 8 else None
    prL = exact_primal(l, d, tau, L) if l <= 8 else None
    # is the pruned set inclusion-minimal? (try deleting single nodes)
    red = sum(1 for e in Bp if exact_dual(l, d, tau, [f for f in Bp if f != e])[1]) if len(Bp) <= 60 else None
    print(f"(l,d)=({l},{d}): bound={bound(l,d)} layered={len(L)} (primal ok {prL}) pruned={len(Bp)} dual_ok={ok} primal_ok={pr} "
          f"per_dir={per_dir(l,Bp)} max={max(per_dir(l,Bp))} layered_max={max(per_dir(l,L))} single-deletable={red}", flush=True)
