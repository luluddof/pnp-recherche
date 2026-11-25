# descending-p scan for fixed (N,q): for p from p0 down, R restarts; 'EXACT?' if tail<1e-18 and gap>1e6 (then gauge-check)
import numpy as np, sys, time, json
from scipy.optimize import minimize
from vp11 import fg, diag, swap
import gauge
N,q,p0,R,seed,tmax=map(int,sys.argv[1:7]); it=int(sys.argv[7]) if len(sys.argv)>7 else 3000
F=swap(N); rng=np.random.default_rng(seed); t0=time.time(); p=p0
while p>0 and time.time()-t0<tmax:
    found=False; bestt=1e9
    for r in range(R):
        if time.time()-t0>tmax: break
        z0=rng.standard_normal(2*q*N*N)*rng.choice([.4,.7,1.])
        z=minimize(lambda z: fg(z,F,N,p,q,'f2',0,0),z0,jac=True,method='L-BFGS-B',options=dict(maxiter=it,ftol=1e-17,gtol=1e-14,maxcor=30)).x
        d=diag(z,N,p,q); bestt=min(bestt,d['tail'])
        if d['tail']<1e-12:
            s=N*N; G2,H2=gauge.run(z[:q*s].reshape(q,N,N),z[q*s:].reshape(q,N,N),N,1)
            zb=np.concatenate([G2.ravel(),H2.ravel()])
            zb=minimize(lambda z: fg(z,F,N,p,q,'f2',0,0),zb,jac=True,method='L-BFGS-B',options=dict(maxiter=4000,ftol=1e-30,gtol=1e-18,maxcor=50)).x
            G2,H2=zb[:q*s].reshape(q,N,N),zb[q*s:].reshape(q,N,N)
            K=np.einsum('kxu,kyv->xyuv',G2,H2).reshape(s,s); sv=np.linalg.svd(F.reshape(s,s)-K,compute_uv=False)
            print(f"N={N} q={q} p={p} r{r} EXACT? tail={d['tail']:.1e} gap={d['gap']:.1e} coef={d['coef']:.1f} Qmin={gauge.Qop(G2,H2):.2f} polished: gap={sv[p-1]/sv[p]:.1e} sv_p+1={sv[p]:.1e} coef={max(abs(G2).max(),abs(H2).max()):.2f} trTh={d['trTh']:.2f} t={time.time()-t0:.0f}",flush=True)
            peff=int(np.sum(sv>1e-7*sv[0])); print(f'   effective numerical rank {peff}, sv[peff-1]/sv[peff]={sv[peff-1]/sv[peff]:.1e}',flush=True)
            np.save(f"ex_{N}_{peff}_{q}_s{seed}.npy",np.concatenate([G2.ravel(),H2.ravel()])); found=True; p=peff; break
    if not found:
        print(f"N={N} q={q} p={p}: not found in {r+1} restarts, best tail {bestt:.2e} t={time.time()-t0:.0f}",flush=True); break
    p-=1
