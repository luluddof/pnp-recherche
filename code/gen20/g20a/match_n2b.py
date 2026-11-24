import itertools, numpy as np, sympy, sys, json
sys.path.insert(0,'/home/user/pnp-recherche/code/gen18/g18h')
from dmatch import matchings, tensor
exec(open('match_n2.py').read().split("rng=np.random")[0].split("ms=list")[0])
ms=list(matchings(list(range(8)),lambda a,b:True)); K=len(ms)
T=np.array([tensor(m,n) for m in ms]).astype(float)
def flat(v,j):
    A=v.reshape([2]*8); rows=[c for c in range(8) if not (c>>j)&1]; cols=[c for c in range(8) if (c>>j)&1]
    return A.transpose([7-c for c in rows]+[7-c for c in cols]).reshape(16,16)
def rk(v,j): return np.linalg.matrix_rank(flat(v,j),tol=1e-8)
mr=np.array([[rk(T[i],j) for j in range(3)] for i in range(K)])
rng=np.random.default_rng(1); Tt=T.T; found={}
target=[(1,2,3),(2,2,2),(1,1,5),(1,1,4),(1,3,3),(1,2,4)]
for trial in range(40000):
    w=rng.random(K)*np.exp(-0.7*mr.min(1)); order=np.argsort(-w); sel=[]
    for i in order:
        if np.linalg.matrix_rank(Tt[:,sel+[i]])>len(sel): sel.append(i)
        if len(sel)==14: break
    c=np.linalg.lstsq(Tt[:,sel],Te.astype(float),rcond=None)[0]
    for rep in range(6):
        asg=[]
        A=[np.zeros(256) for _ in range(3)]
        for ci,i in zip(c,sel):
            r=mr[i]; cand=np.flatnonzero(r==r.min()) if rep==0 else np.flatnonzero(r<=np.sort(r)[1])
            j=int(rng.choice(cand)); A[j]+=ci*T[i]; asg.append(j)
        rr=[rk(A[j],j) for j in range(3)]
        key=tuple(sorted(rr))
        if key in target and key not in found:
            found[key]=1
            # exact verification
            M=sympy.Matrix(T[sel].T.astype(int).tolist()); b=sympy.Matrix(Te.tolist())
            cs=M.solve_least_squares(b) if False else M.gauss_jordan_solve(b)[0]
            assert M*cs==b
            Aex=[sympy.zeros(256,1) for _ in range(3)]
            for k,(i,j) in enumerate(zip(sel,asg)): Aex[j]+=cs[k]*sympy.Matrix(T[i].astype(int).tolist())
            tot=Aex[0]+Aex[1]+Aex[2]; assert tot==b
            ex=[]
            for j in range(3):
                v=np.array([sympy.nsimplify(x) for x in Aex[j]],dtype=object)
                F=sympy.Matrix(flat(v,j).tolist()); ex.append(F.rank())
            print('trial',trial,'float',rr,'EXACT ranks',ex,'coeffs',sorted(set(cs)),'matchings',[(ms[i],int(j)) for i,j in zip(sel,asg)],flush=True)
            json.dump({'sel':[int(i) for i in sel],'asg':asg,'coef':[str(x) for x in cs],'ranks':ex},open(f'exact_n2_{"_".join(map(str,ex))}.json','w'))
    if len(found)==len(target): break
print('found',found)
