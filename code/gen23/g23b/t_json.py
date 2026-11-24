import glob, sys, json, time
from core import *
for f in sorted(glob.glob('/home/user/pnp-recherche/code/gen22/g22b/shape_*.json')):
    N, col, pal, sg = load_json(f)
    M = classes(N, col)
    nS = sum(1 for k in M if pal[k]=='S'); nT = sum(1 for k in M if pal[k]=='T'); nU = sum(1 for k in M if pal[k]=='U')
    t=time.time(); ok, r, c = ranks(N, col, pal, sg)
    print(f.split('/')[-1], 'N',N,'S,T,U',nS,nT,nU,'tensor id',identity_tensor_ok(N,col,pal,sg),'pred',predicate(N,col,pal,sg),
          'R+C=X',ok,'rk R',r,'rk C',c, 'lemma', key_lemma_check(N,col,pal), 'sizes', sorted(len(v) for v in M.values()),
          'bound slack', nS*gstar(nT,N)+nT*gstar(nS,N)-N*N if nU==0 else None, f'{time.time()-t:.1f}s')
