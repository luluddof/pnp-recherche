# Independent end-to-end test of the mod-3 evaluator (E2 + per-point block design) against brute force.
# Pipeline: labels -> copies -> rows P_{f,v}(T) over GF(3^k) -> per-node measurement matrices mod 3 ->
#           separated evaluation (Lemma 3.2 mod 3) -> Taylor data via Lemma 3.3 -> Hermite/CRT recovery -> H(0).
import numpy as np, itertools, math, random, sys, time
from gf import GF
from design import Design, eigvecs

def gen_formula(n, ncl, rng, maxocc=3):
    while True:
        occ=[0]*n; cls=[]; ok=True
        for _ in range(ncl):
            cand=[v for v in range(n) if occ[v]<maxocc]
            if len(cand)<3: ok=False; break
            vs=rng.sample(cand,3)
            for v in vs: occ[v]+=1
            cls.append([(v,rng.randrange(2)) for v in vs])   # (var, neg?)
        if ok: return cls

class VF:
    """vectorised GF ops"""
    def __init__(s,F):
        s.F=F; s.ADD=F.ADD; s.MUL=F.MUL; s.NEG=F.NEG; s.k=F.k
        s.DIG=np.array([[(a//3**i)%3 for i in range(F.k)] for a in range(F.q)],dtype=np.int64)
        s.POW=np.array([3**i for i in range(F.k)],dtype=np.int64)
    def add(s,a,b): return s.ADD[a,b]
    def mul(s,a,b): return s.MUL[a,b]
    def sub(s,a,b): return s.ADD[a,s.NEG[b]]
    def lincomb_int(s, Mint, A, axis):
        """integer matrix Mint (rows x k) applied along 'axis' of GF array A (sum with small integer coefficients)"""
        Dg=s.DIG[A]                                   # (..., kdig)
        Dg=np.moveaxis(Dg,axis,-2)                    # (..., k, kdig)
        R=np.einsum('rk,...kd->...rd',Mint,Dg)%3      # (..., rows, kdig)
        R=R@s.POW
        return np.moveaxis(R,-1,axis)
    def pmul(s,A,B,J):
        La=A.shape[-1]; Lb=B.shape[-1]; L=min(J,La+Lb-1)
        shp=np.broadcast_shapes(A.shape[:-1],B.shape[:-1])
        out=np.zeros(shp+(L,),dtype=np.int64)
        for j in range(Lb):
            w=min(La,L-j)
            if w<=0: break
            out[...,j:j+w]=s.ADD[out[...,j:j+w], s.MUL[A[...,:w],B[...,j:j+1]]]
        return out

def taylor_shift(F,p,tau):
    """coefficients of p(tau+z) in z (Hasse-Taylor): sum_k p_k C(k,j) tau^(k-j)"""
    out=[]
    for j in range(len(p)):
        acc=0
        for k in range(j,len(p)):
            c=math.comb(k,j)%3
            if c==0 or p[k]==0: continue
            term=F.mul(c,p[k]); pw=1
            for _ in range(k-j): pw=F.mul(pw,tau)
            acc=F.add(acc,F.mul(term,pw))
        out.append(acc)
    return out

def solve(F,A,b):
    """Gaussian elimination over GF; returns x or None if singular"""
    n=len(A); A=[row[:]+[b[i]] for i,row in enumerate(A)]
    for c in range(n):
        piv=next((i for i in range(c,n) if A[i][c]),None)
        if piv is None: return None
        A[c],A[piv]=A[piv],A[c]
        iv=F.inv(A[c][c]); A[c]=[F.mul(iv,x) for x in A[c]]
        for i in range(n):
            if i!=c and A[i][c]:
                f=A[i][c]; A[i]=[F.sub(x,F.mul(f,y)) for x,y in zip(A[i],A[c])]
    return [A[i][n] for i in range(n)]

def rowbasis_mod3(rows, ncols):
    """row basis (reduced echelon) mod 3 of a list of int vectors"""
    basis=[]; piv=[]
    for v in rows:
        v=np.array(v,dtype=np.int64)%3
        for b,p in zip(basis,piv):
            if v[p]: v=(v-v[p]*b)%3
        nz=np.nonzero(v)[0]
        if len(nz)==0: continue
        p=nz[0]
        if v[p]==2: v=(2*v)%3
        for i in range(len(basis)):
            if basis[i][p]: basis[i]=(basis[i]-basis[i][p]*v)%3
        basis.append(v); piv.append(p)
    return np.array(basis,dtype=np.int64).reshape(len(basis),ncols)

def run(n, ncl, r, Q, m, shift, k=2, seed=0, pool=None, weights=True, verbose=True, maxpairs=None):
    rng=random.Random(seed); F=GF(k); V=VF(F); t0=time.time()
    cls=gen_formula(n,ncl,rng)
    taus=rng.sample(range(1,F.q),Q)
    des=Design(F,Q,m,shift,taus)
    if pool:
        pa=[tuple(rng.randrange(r) for _ in range(Q)) for _ in range(pool)]
        pd=[tuple(rng.randrange(r) for _ in range(Q)) for _ in range(pool)]
        a=[rng.choice(pa) for _ in range(n)]; d=[rng.choice(pd) for _ in range(ncl)]
    else:
        a=[tuple(rng.randrange(r) for _ in range(Q)) for _ in range(n)]
        d=[tuple(rng.randrange(r) for _ in range(Q)) for _ in range(ncl)]
    w=[rng.choice([e for e in range(1,F.q) if e!=2]) if weights else 1 for _ in range(n)]
    pairs=[(f,v) for f,c in enumerate(cls) for (v,_) in c]
    label={(f,v):tuple((a[v][s],d[f][s]) for s in range(Q)) for (f,v) in pairs}
    cnt={}; copy={}
    for pr in pairs:
        copy[pr]=cnt.get(label[pr],0); cnt[label[pr]]=copy[pr]+1
    h=max(cnt.values())
    zs=list(itertools.product((0,1),repeat=Q))
    def pi(zs_,p,p2):   # entry of (I + (-1)^z Sw)/2 mod 3 at ((a,d),(a',d'))
        val=(1 if p==p2 else 0) + ((-1)**zs_)*(1 if p==(p2[1],p2[0]) else 0)
        return (2*val)%3
    def Pentry(u,y):
        poly=[0]*m
        for z in zs:
            c=1
            for s in range(Q):
                c=(c*pi(z[s],u[s],y[s]))%3
                if c==0: break
            if c: poly=F.padd(poly, F.pscal(c, des.Pz(z)+[0]*(m-len(des.Pz(z)))))
        return poly[:m]
    rows={}
    for pr in pairs:
        R=[[0]*m for _ in range(n)]
        for pr2 in pairs:
            if copy[pr2]==copy[pr]:
                e=Pentry(label[pr],label[pr2])
                R[pr2[1]]=F.padd(R[pr2[1]],e)
        rows[pr]=R
        for v2 in range(n):
            assert R[v2][0]==(1 if v2==pr[1] else 0), "P_{f,v}(0) != delta_v"
    degrow=max(F.pdeg(R[v2]) for R in rows.values() for v2 in range(n))
    dH=3*ncl*degrow if degrow>0 else 0
    q=dH//Q; J=q+1
    assert Q*J>dH
    # ---------- brute force ----------
    X=np.array(list(itertools.product((0,1),repeat=n)),dtype=np.int64)   # (2^n, n)
    sat=np.ones(len(X),dtype=bool)
    for c in cls:
        fals=np.ones(len(X),dtype=bool)
        for (v,neg) in c: fals&=(X[:,v]==(1 if neg else 0))
        sat&=~fals
    nsat=int(sat.sum())
    wx=np.ones(len(X),dtype=np.int64)
    for i in range(n): wx=V.mul(wx, np.where(X[:,i]==1,w[i],1))
    truth=int(V.lincomb_int(sat.astype(np.int64)[None,:], wx, 0)[0])    # sum_{x sat} w^x in GF
    ev=eigvecs(r)
    def psi_val(choice,u):
        val=1
        for s in range(Q):
            val*=int(ev[choice[s]][0][u[s][0]*r+u[s][1]])
            if val==0: return 0
        return val
    def brute_node(t):
        """H(tau_t+z) mod z^J by brute force over all x"""
        prod=np.zeros((len(X),J),dtype=np.int64); prod[:,0]=wx
        for f,c in enumerate(cls):
            fal=np.zeros((len(X),1),dtype=np.int64); fal[:,0]=1
            for (v,neg) in c:
                sh=np.array([taylor_shift(F,rows[(f,v)][i],taus[t]) for i in range(n)],dtype=np.int64)  # (n,m)
                lin=np.stack([V.lincomb_int(X, sh[:,j], 0) for j in range(m)],axis=-1)   # (2^n, m)
                if not neg:
                    lin=V.NEG[lin]; lin[:,0]=V.add(lin[:,0],1)
                fal=V.pmul(fal,lin,J)
            phi=V.NEG[fal]; phi[:,0]=V.add(phi[:,0],1)
            prod=V.pmul(prod,phi,J)
        return V.lincomb_int(np.ones((1,len(X)),dtype=np.int64), prod, 0)[0]
    stats=dict(n=n,ncl=ncl,r=r,Q=Q,m=m,shift=shift,h=h,M=len(pairs),degrow=degrow,dH=dH,q=q,J=J,nsat=nsat)
    Elist=[E for s in range(min(q,n)+1) for E in itertools.combinations(range(n),s)]
    npairs=sum(2**len(E) for E in Elist); stats['pairs_per_node']=npairs
    data=[]; kernel_nonvac=0; kernel_checks=0; tuples_info=[]
    for t in range(Q):
        S=[[v for v in range(n) if a[v][t]==g] for g in range(r)]
        # measurement matrices
        bad=[ch for ch in itertools.product(range(len(ev)),repeat=Q)
             if sum(ev[c][1] for s,c in enumerate(ch) if s!=t) in des.tops]
        mrows=[]
        for cidx in range(h):
            occ=[pr for pr in pairs if copy[pr]==cidx]
            for ch in bad:
                v_=np.zeros(n,dtype=np.int64)
                for pr in occ: v_[pr[1]]+=psi_val(ch,label[pr])
                if v_.any(): mrows.append(v_)
        U=[]; clsid=[]; ncls=[]; rks=[]
        for g in range(r):
            Bg=rowbasis_mod3([row[S[g]] for row in mrows], len(S[g]))
            Ug=np.array(list(itertools.product((0,1),repeat=len(S[g]))),dtype=np.int64).reshape(2**len(S[g]),len(S[g]))
            meas=(Ug@Bg.T)%3 if len(Bg) else np.zeros((len(Ug),0),dtype=np.int64)
            keys={}; ids=np.zeros(len(Ug),dtype=np.int64)
            for i,row in enumerate(map(tuple,meas)): ids[i]=keys.setdefault(row,len(keys))
            U.append(Ug); clsid.append(ids); ncls.append(len(keys)); rks.append(len(Bg))
        tuples_info.append((rks,[len(s_) for s_ in S],ncls))
        # per factor/variable tables
        sig=[d[f][t] for f in range(ncl)]
        sh={pr:np.array([taylor_shift(F,rows[pr][i],taus[t]) for i in range(n)],dtype=np.int64) for pr in pairs} # (n,m)
        tabs={}
        for pr in pairs:
            g=sig[pr[0]]; rho=sh[pr][:,0]
            tb=[]
            for h_ in range(r):
                full=V.lincomb_int(U[h_], rho[S[h_]], 0) if len(S[h_]) else np.zeros(1,dtype=np.int64)
                if h_==g: tb.append(full)
                else:
                    cv=np.full(ncls[h_],-1,dtype=np.int64)
                    for i,cid in enumerate(clsid[h_]):
                        if cv[cid]<0: cv[cid]=full[i]
                        assert cv[cid]==full[i], ("KERNEL CONDITION VIOLATED", t, pr, h_)
                    kernel_checks+=1
                    if len(U[h_])>ncls[h_] and len(set(full.tolist()))>1: kernel_nonvac+=1
                    tb.append(cv)
            tabs[pr]=tb
        wtab=[]
        for g in range(r):
            wg=np.ones(len(U[g]),dtype=np.int64)
            for j,i in enumerate(S[g]): wg=V.mul(wg,np.where(U[g][:,j]==1,w[i],1))
            wtab.append(wg)
        onehot=[np.array([[1 if clsid[g][i]==c else 0 for i in range(len(U[g]))] for c in range(ncls[g])],dtype=np.int64) for g in range(r)]
        def shape_for(g,arr,axis_h):   # reshape 1-d arr to broadcast along axis_h among r axes
            shp=[1]*r; shp[axis_h]=len(arr); return arr.reshape(shp)
        Ht=np.zeros(J,dtype=np.int64)
        for E in Elist:
            Eset=set(E)
            for su in range(len(E)+1):
                for Usub in itertools.combinations(E,su):
                    sign=(-1)**(len(E)-su)
                    Ssum=[]
                    for g in range(r):
                        dims=[ncls[h_] for h_ in range(r)]; dims[g]=len(U[g])
                        Ag=np.zeros(tuple(dims)+(1,),dtype=np.int64); Ag[...,0]=shape_for(g,wtab[g],g)
                        mask=np.ones(len(U[g]),dtype=np.int64)
                        for j,i in enumerate(S[g]):
                            if i in Eset: mask*=U[g][:,j]
                        for f in range(ncl):
                            if sig[f]!=g: continue
                            fal=np.zeros(tuple(dims)+(1,),dtype=np.int64); fal[...,0]=1
                            for (v,neg) in cls[f]:
                                tb=tabs[(f,v)]
                                base=np.zeros(tuple(dims),dtype=np.int64)
                                for h_ in range(r): base=V.add(base, shape_for(g,tb[h_],h_))
                                lin=np.zeros(tuple(dims)+(m,),dtype=np.int64); lin[...,0]=base
                                if Usub:
                                    corr=V.lincomb_int(np.ones((1,len(Usub)),dtype=np.int64), sh[(f,v)][list(Usub),:], 0)[0]  # (m,)
                                    lin[...,1:]=corr[1:]
                                if not neg:
                                    lin=V.NEG[lin]; lin[...,0]=V.add(lin[...,0],1)
                                fal=V.pmul(fal,lin,J)
                            phi=V.NEG[fal]; phi[...,0]=V.add(phi[...,0],1)
                            Ag=V.pmul(Ag,phi,J)
                        if Ag.shape[-1]<J:
                            Ag=np.concatenate([Ag,np.zeros(Ag.shape[:-1]+(J-Ag.shape[-1],),dtype=np.int64)],axis=-1)
                        Mg=onehot[g]*mask[None,:]
                        Ssum.append(V.lincomb_int(Mg, Ag, g))      # axis g now indexes own class
                    tot=Ssum[0]
                    for g in range(1,r): tot=V.pmul(tot,Ssum[g],J)
                    flat=tot.reshape(-1,tot.shape[-1])
                    val=V.lincomb_int(np.ones((1,len(flat)),dtype=np.int64), flat, 0)[0]
                    if val.shape[0]<J: val=np.concatenate([val,np.zeros(J-val.shape[0],dtype=np.int64)])
                    Ht=V.add(Ht, val if sign==1 else V.NEG[val])
        bf=brute_node(t)
        assert (Ht==bf).all(), ("Taylor data mismatch at node", t, Ht, bf)
        data.append(Ht)
    # ---------- recovery (confluent Vandermonde with Hasse derivatives) ----------
    N=Q*J; A=[]; b=[]
    for t in range(Q):
        pw=[1]
        for _ in range(N): pw.append(F.mul(pw[-1],taus[t]))
        for j in range(J):
            A.append([F.mul(math.comb(k_,j)%3, pw[k_-j]) if k_>=j else 0 for k_ in range(N)]); b.append(int(data[t][j]))
    coeffs=solve(F,A,b)
    assert coeffs is not None, "Hermite matrix singular over GF"
    assert all(c==0 for c in coeffs[dH+1:]), "recovered degree exceeds d_H"
    got=coeffs[0]
    stats.update(truth=truth, got=got, nsat_mod3=nsat%3, ok=(got==truth), kernel_checks=kernel_checks,
                 kernel_nonvac=kernel_nonvac, tuples=tuples_info, time=round(time.time()-t0,1),
                 degH=max([i for i,c in enumerate(coeffs) if c],default=-1))
    if not weights: assert truth==nsat%3
    assert got==truth, "FINAL MISMATCH"
    if verbose: print(stats, flush=True)
    return stats

if __name__=="__main__":
    mode=sys.argv[1]
    if mode=="quick":
        run(6,2,2,3,3,3,seed=1); run(6,2,2,3,3,3,seed=2,weights=False)
    else:
        cfgs=eval(sys.argv[2])
        for c in cfgs: run(**c)
