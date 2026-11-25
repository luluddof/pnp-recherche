# D = inf (full Hom = B6 Test A in operator form) for larger h: per sector, per W: quotient by Rho(W), then add kappa blocks
# in increasing W' until all rhs lie in the span.  Exact mod P ranks.  Also reports the deficiency (#rhs missing) per W'.
# usage: python3 hinf.py h Wlist Wpmax [P]
import sys, time
sys.dont_write_bytecode = True
import numpy as np
from sscan import Setup, rowbasis, rank, fl
from fcomp import mm

def reduce_mod(V, Rr, piv, nonpiv, P):
    """V (k, m) modulo row space of rref Rr (pivots piv) -> coordinates on nonpivot columns"""
    if Rr.shape[0] == 0: return V[:, nonpiv] % P
    X = (V[:, nonpiv] - mm(V[:, piv] % P, Rr[:, nonpiv], P)) % P
    return X

def main(h, Ws, Wpmax, P=16777213):
    t0 = time.time(); Wmax = max(Ws); st = Setup(h, Wmax, Wpmax, P); n = st.n
    print(f"h={h} n={n} Ws={Ws} W'max={Wpmax}: setup {time.time()-t0:.0f}s", flush=True)
    full = {}
    for side, rs in (('rho', st.rr), ('kap', st.rk)):
        for r in rs:
            for pi in (0, 1):
                m = st.hmask[(side, r)]; idx = np.nonzero(m == pi)[0]
                E = np.zeros((len(idx), len(m)), dtype=np.int64); E[np.arange(len(idx)), idx] = 1; full[(side, r, pi)] = E
    res = {W: [] for W in Ws}
    for pi in (0, 1):
        R = st.rhs(pi); m = len(st.rows[pi])
        ck = {}
        for r in st.rk:
            ck[r] = rowbasis(st.cols('kap', r, full[('kap', r, pi)], pi), P)
        print(f" sector {pi}: rows {m}, kappa blocks {[ck[r].shape[0] for r in st.rk]}  {time.time()-t0:.0f}s", flush=True)
        Rho = np.zeros((0, m), dtype=np.int64); Wdone = -1
        for W in sorted(Ws):
            for Wr in range(Wdone + 1, W + 1):
                r = Wr - h
                if r in st.rr:
                    blk = rowbasis(st.cols('rho', r, full[('rho', r, pi)], pi), P)
                    Rho = rowbasis(np.concatenate([Rho, blk]), P)
            Wdone = W
            piv = np.array([int(np.nonzero(row)[0][0]) for row in Rho], dtype=np.int64)
            nonpiv = np.setdiff1d(np.arange(m), piv)
            Rred = rowbasis(reduce_mod(R, Rho, piv, nonpiv, P), P)
            Kb = np.zeros((0, len(nonpiv)), dtype=np.int64); found = None; prof = []
            for Wp in range(-1, Wpmax + 1):
                if Wp >= 0:
                    r = h - Wp
                    if r in ck and ck[r].shape[0]:
                        Kb = rowbasis(np.concatenate([Kb, reduce_mod(ck[r], Rho, piv, nonpiv, P)]), P)
                tot = rank(np.concatenate([Kb, Rred]), P) if Rred.shape[0] else Kb.shape[0]
                defi = tot - Kb.shape[0]; prof.append(defi)
                if defi == 0: found = Wp; break
                if 4 * h - 1 - W - Wp < 0: break
            res[W].append(found)
            print(f"  sector {pi} W={W}: rank Rho={Rho.shape[0]} quotient dim {len(nonpiv)}; min W'={found}; deficiency profile from W'=-1: {prof}  {time.time()-t0:.0f}s", flush=True)
    for W in Ws:
        l = res[W]
        if None in l: print(f"h={h} D=inf W={W}: min W' > {Wpmax}")
        else:
            Wp = max(l); print(f"h={h} D=inf W={W}: min W'={Wp}  k={4*h-1-W-Wp}")

if __name__ == "__main__":
    h = int(sys.argv[1]); Ws = [int(a) for a in sys.argv[2].split(',')]; Wpmax = int(sys.argv[3])
    P = int(sys.argv[4]) if len(sys.argv) > 4 else 16777213
    main(h, Ws, Wpmax, P)
