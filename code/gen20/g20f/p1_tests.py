import numpy as np, itertools, random
from p1_lib import *
rng=np.random.default_rng(1); random.seed(1)
# (1a) cut-conditioning identity: exact tensor check
def check_cc(wires,n,lab):
    M=len(wires); bad=0
    cut=[e for e,(h,t) in enumerate(wires) if lab[h]!=lab[t]]
    for bits in itertools.product((0,1),repeat=2*M):
        Y=bits[:M]; X=bits[M:]; I=int(all(Y[e]==X[e] for e in range(M)))
        s=0
        for z in itertools.product((0,1),repeat=len(cut)):
            zz=dict(zip(cut,z)); prod=1
            for g in set(lab):
                f=1
                for e,(h,t) in enumerate(wires):
                    if e in zz:
                        if lab[h]==g and Y[e]!=zz[e]: f=0      # row end Y_e at h
                        if lab[t]==g and X[e]!=zz[e]: f=0      # column end X_e at t
                    elif lab[h]==g and Y[e]!=X[e]: f=0
                prod*=f
            s+=prod
        bad+= (s!=I)
    return bad
tot=0
for trial in range(30):
    n=random.randint(2,5); M=random.randint(1,4); r=random.randint(1,4)
    wires=[(random.randrange(n),random.randrange(n)) for _ in range(M)]
    lab=np.array([random.randrange(r) for _ in range(n)])
    tot+=check_cc(wires,n,lab)
print('(1a) cut-conditioning identity: mismatching entries over 30 random instances (incl. loops, empty groups):',tot)
# (1c) Lemma Q / Theorem E, brute force A''
viol_r2=viol_gen=0; cnt=0; tight=[]
for trial in range(600):
    n=random.randint(2,7); M=random.randint(1,9); u=random.randint(1,4); r=random.choice([2,2,3,4,5])
    wires=[(random.randrange(n),random.randrange(n)) for _ in range(M)]
    # adversarial: skewed group sizes, empty groups allowed
    parts=[np.array([min(r-1,int(rng.exponential(0.7))) for _ in range(n)]) for _ in range(u)]
    A,B,X=Aexact(wires,n,parts); mi=Mi(wires,parts); cnt+=1
    if r==2 and A < min(mi)/u-u-1e-9: viol_r2+=1
    if A < min(mi)/(2*u)-u-1e-9: viol_gen+=1
    tight.append((A - (min(mi)/u if r==2 else min(mi)/(2*u)), r))
print('(1c) instances',cnt,' violations r=2 bound:',viol_r2,' general bound:',viol_gen)
print('     min of A - minMi/u (r=2):',min(t for t,r in tight if r==2),' min of A - minMi/(2u) (r>=3):',min(t for t,r in tight if r>2))
