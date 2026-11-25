import sys, json, numpy as np, mc
def load(fn, kE, kC):
    d = json.load(open(fn)); E = [np.array(m, dtype=float).astype(int) for m in d[kE]]; C = [np.array(m, dtype=float).astype(int) for m in d[kC]]
    return E, C
for fn, kE, kC in [(a.split(':')[0], a.split(':')[1], a.split(':')[2]) for a in sys.argv[1:]]:
    E, C = load(fn, kE, kC); N = E[0].shape[0]; q = len(E)
    K = mc.Kmat(N, E, C); S = mc.swap(N); Th = K @ S
    n = N * N
    # peeling: core
    alive = np.ones(n, bool)
    A = (Th != 0)
    while True:
        indeg = (A[np.ix_(alive, alive)]).sum(1); outdeg = (A[np.ix_(alive, alive)]).sum(0)
        idx = np.nonzero(alive)[0]; kill = idx[(indeg == 0) | (outdeg == 0)]
        if len(kill) == 0: break
        alive[kill] = False
    core = np.nonzero(alive)[0]
    import flint
    rk = flint.fmpz_mat((np.eye(n, dtype=int) - Th).tolist()).rank()
    rkc = flint.fmpz_mat((np.eye(len(core), dtype=int) - Th[np.ix_(core, core)]).tolist()).rank() if len(core) else 0
    diag = np.diag(Th)
    print(f"{fn.split('/')[-1]}: N={N} q={q} rank={rk} |core|={len(core)} rank(I-Th_core)={rkc}  diag(Th): +1:{(diag==1).sum()} -1:{(diag==-1).sum()} 0:{(diag==0).sum()}  nnz(Th)={int((Th!=0).sum())} ranks E {[int(abs(m).sum()) for m in E]} C {[int(abs(m).sum()) for m in C]}")
