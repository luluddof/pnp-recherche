# Faster joint-greedy search for the G-invariant model (see wreath.py for model).
import numpy as np, sys, time, json
from wreath import build, mat, coeffs, center_idempotents, rk
def setup(ns, seed=0):
    rng = np.random.default_rng(seed); t0 = time.time()
    S = build(ns); N, r = S['N'], S['r']; D = N*N
    E0, zc = center_idempotents(S, rng)
    SW = np.zeros((D, D)); a = np.arange(D); SW[(a % N) * N + a // N, a] = 1
    E = []
    for e in E0:
        for sg in (1, -1):
            f = e @ (np.eye(D) + sg * SW) / 2
            if np.trace(f) > 0.5: E.append(f)
    w = [int(round(np.trace(e))) for e in E]
    B = [mat(np.eye(r)[o], S) for o in range(r)]
    LR = [np.array([coeffs(Bo @ e, S) for Bo in B]).T for e in E]
    sR = [coeffs(SW @ e, S) for e in E]
    P = np.zeros((r, r)); P[S['piR'], np.arange(r)] = 1
    LC = [L @ P for L in LR]
    print(f"ns={ns} N={N} |orbits|={r} center dim={zc} pieces={len(E)} weights={w} ({time.time()-t0:.1f}s)", flush=True)
    return S, SW, E, w, LR, sR, LC, rng
def search(ns, NS=2000, seed=0, out=None):
    S, SW, E, w, LR, sR, LC, rng = setup(ns, seed)
    N, r = S['N'], S['r']; D = N*N; nI = len(E); t0 = time.time()
    best = {}
    cons = [('R', i) for i in range(nI)] + [('C', i) for i in range(nI)]
    for it in range(NS):
        lam = rng.random()
        keys = np.array([(rng.random() * (lam if t == 'R' else 1 - lam)) * (1 + w[i]) ** rng.uniform(0, 2) for t, i in cons])
        order = np.argsort(-keys)
        c0 = np.zeros(r); Nb = np.eye(r); used = []
        for j in order:
            t, i = cons[j]
            L = LR[i] if t == 'R' else LC[i]; s = sR[i] if t == 'R' else np.zeros(r)
            M = L @ Nb; rhs = s - L @ c0
            if Nb.shape[1] == 0:
                if np.linalg.norm(rhs) < 1e-7: used.append(j)
                continue
            u, sv, vt = np.linalg.svd(M, full_matrices=True)
            rnk = int(np.sum(sv > 1e-9 * max(1.0, sv[0] if len(sv) else 1)))
            z = vt[:rnk].T @ ((u[:, :rnk].T @ rhs) / sv[:rnk]) if rnk else np.zeros(Nb.shape[1])
            if np.linalg.norm(M @ z - rhs) > 1e-7 * (1 + np.linalg.norm(rhs)): continue
            c0 = c0 + Nb @ z; Nb = Nb @ vt[rnk:].T; used.append(j)
        c = c0 + (Nb @ rng.standard_normal(Nb.shape[1]) if Nb.shape[1] else 0)
        K = mat(c, S)
        p = rk(SW - K); q = rk(K.reshape(N, N, N, N).transpose(0, 2, 1, 3).reshape(D, D))
        if q not in best or p < best[q][0]:
            best[q] = (p, [cons[j] for j in used], c.tolist())
    fr = []; bp = 10**9
    for q in sorted(best):
        if best[q][0] < bp: bp = best[q][0]; fr.append((best[q][0], q))
    print("frontier (p,q):", fr, f"({time.time()-t0:.1f}s)", flush=True)
    if out: json.dump({str(q): [best[q][0], best[q][1], best[q][2]] for q in best}, open(out, 'w'))
    return fr
if __name__ == '__main__':
    ns = list(map(int, sys.argv[1].split(','))); NS = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    search(ns, NS, out=f"wreath2_{'_'.join(map(str,ns))}.json")
