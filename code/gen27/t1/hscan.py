# Large-h version of sscan: rho mode spans S_r(D) (sources V_h only), kappa spans by adjoint (validated), and per W an
# incremental quotient scan over W' (memory O(rows^2) instead of storing all prefix spans).  Exact mod P, no projection.
# usage: python3 hscan.py h Dlist Ws Wmax [P]     (Dlist may contain 'inf')
import sys, time, pickle, os
sys.dont_write_bytecode = True
import numpy as np
from sscan import Setup, rowbasis, rank, adjoint_fill
from hinf import reduce_mod
from appb_exact import partitions

def scan_inc(st, S, Ws, label, log, kmin=None):
    h, P = st.h, st.P; t0 = time.time(); res = {W: [] for W in Ws}; prof_all = {}
    for pi in (0, 1):
        R = st.rhs(pi); m = len(st.rows[pi])
        ck = {r: rowbasis(st.cols('kap', r, S[('kap', r, pi)], pi), P) for r in st.rk}
        Rho = np.zeros((0, m), dtype=np.int64); Wdone = -1
        for W in sorted(Ws):
            for Wr in range(Wdone + 1, W + 1):
                r = Wr - h
                if r in st.rr:
                    Rho = rowbasis(np.concatenate([Rho, st.cols('rho', r, S[('rho', r, pi)], pi)]), P)
            Wdone = W
            piv = np.array([int(np.nonzero(row)[0][0]) for row in Rho], dtype=np.int64)
            nonpiv = np.setdiff1d(np.arange(m), piv)
            Rred = rowbasis(reduce_mod(R, Rho, piv, nonpiv, P), P)
            Kb = np.zeros((0, len(nonpiv)), dtype=np.int64); found = None; prof = []
            for Wp in range(-1, st.Wpmax + 1):
                if Wp >= 0:
                    r = h - Wp
                    if r in ck and ck[r].shape[0]:
                        Kb = rowbasis(np.concatenate([Kb, reduce_mod(ck[r], Rho, piv, nonpiv, P)]), P)
                if kmin is not None and 4 * h - 1 - W - Wp > kmin + 6 and Wp < st.Wpmax: prof.append('.'); continue
                defi = (rank(np.concatenate([Kb, Rred]), P) - Kb.shape[0]) if Rred.shape[0] else 0
                prof.append(defi)
                if defi == 0: found = Wp; break
            res[W].append(found); prof_all[(pi, W)] = prof
            print(f"   [{label}] sector {pi} W={W}: rank Rho={Rho.shape[0]}/{m}; min W'={found}; deficiency from W'=-1: {prof}  {time.time()-t0:.0f}s", flush=True)
    out = {}
    for W in Ws:
        l = res[W]; out[W] = None if None in l else (max(l), 4 * h - 1 - W - max(l))
    line = f"h={h} D={label}: " + " ".join(f"W{W}:{'-' if out[W] is None else str(out[W][0])+'/k'+str(out[W][1])}" for W in Ws)
    print(line, flush=True); log.write(line + "\n"); log.flush()
    return out, prof_all

def main(h, Dlist, Ws, Wmax, P=16777213):
    t0 = time.time(); st = Setup(h, Wmax, Wmax, P); n = st.n
    log = open(f"hscan_h{h}_P{P}.txt", "a")
    print(f"h={h} n={n} Wmax={Wmax} rows/sector {len(st.rows[0])},{len(st.rows[1])}: setup {time.time()-t0:.0f}s", flush=True)
    keys = [('rho', r, pi) for r in st.rr for pi in (0, 1)]
    fulldim = lambda key: int((st.hmask[key[:2]] == key[2]).sum())
    S = {key: np.zeros((0, len(st.hmask[key[:2]])), dtype=np.int64) for key in keys}
    ck = f"Hck_h{h}_P{P}_W{Wmax}.pkl"; w0 = 0
    if os.path.exists(ck): S, w0 = pickle.load(open(ck, 'rb')); w0 += 1; print(f" resumed, weights < {w0} done", flush=True)
    finite = [D for D in Dlist if D != 'inf']; results = {}
    for w in range(w0, (max(finite) if finite else -1) + 1):
        new = {key: [] for key in keys}
        for mu in partitions(w):
            pi = len(mu) % 2
            for r in st.rr:
                if S[('rho', r, pi)].shape[0] < fulldim(('rho', r, pi)): new[('rho', r, pi)].append(st.O.mat(mu, r, h).ravel())
        for key in keys:
            if new[key]: S[key] = rowbasis(np.concatenate([S[key], np.array(new[key], dtype=np.int64) % P]), P)
        pickle.dump((S, w), open(ck + '.tmp', 'wb')); os.replace(ck + '.tmp', ck)
        msg = f" weight {w}: rho mode-span dim {sum(S[k].shape[0] for k in keys)} / {sum(fulldim(k) for k in keys)}; saturated {sum(S[k].shape[0]==fulldim(k) for k in keys)}/{len(keys)}  {time.time()-t0:.0f}s"
        print(msg, flush=True); log.write(msg + "\n"); log.flush()
        if w in finite:
            S2 = dict(S); adjoint_fill(st, S2, P); results[w] = scan_inc(st, S2, Ws, str(w), log)
    if 'inf' in Dlist:
        Sf = {}
        for side, rs in (('rho', st.rr), ('kap', st.rk)):
            for r in rs:
                for pi in (0, 1):
                    m = st.hmask[(side, r)]; idx = np.nonzero(m == pi)[0]
                    E = np.zeros((len(idx), len(m)), dtype=np.int64); E[np.arange(len(idx)), idx] = 1; Sf[(side, r, pi)] = E
        results['inf'] = scan_inc(st, Sf, Ws, 'inf', log)
    pickle.dump(results, open(f"hscan_h{h}_P{P}_res_{'_'.join(map(str, Dlist))}.pkl", "wb"))
    print(f"total {time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    h = int(sys.argv[1]); Dlist = [a if a == 'inf' else int(a) for a in sys.argv[2].split(',')]
    Ws = [int(a) for a in sys.argv[3].split(',')]; Wmax = int(sys.argv[4]); P = int(sys.argv[5]) if len(sys.argv) > 5 else 16777213
    main(h, Dlist, Ws, Wmax, P)
