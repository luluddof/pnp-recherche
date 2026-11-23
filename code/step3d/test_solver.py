import random, itertools, subprocess
H='/tmp/claude-0/-home-user/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/step3d/howell'
rng=random.Random(1); bad=0; cnt={0:0,1:0}
for t in range(3000):
    q=rng.choice([2,3,4]); nr=rng.randint(1,5); nc=rng.randint(1,5)
    A=[[rng.choice([0,0,1,2,3][:q+2]) % q for _ in range(nc)] for _ in range(nr)]
    if q==4 and rng.random()<.5: A=[[ (2*a)%4 if rng.random()<.6 else a for a in r] for r in A]
    tgt=nc-1
    trip=[f"{i} {j} {A[i][j]}" for i in range(nr) for j in range(nc) if A[i][j]]
    out=subprocess.run([H],input=f"{q} {nr} {nc} {tgt} {len(trip)}\n"+"\n".join(trip)+"\n",capture_output=True,text=True).stdout.strip()=='1'
    e=[0]*nc; e[tgt]=1
    bf=any(all(sum(c[i]*A[i][j] for i in range(nr))%q==e[j] for j in range(nc)) for c in itertools.product(range(q),repeat=nr))
    cnt[bf]+=1
    if out!=bf: bad+=1
print('mismatches',bad,'truth counts',cnt)
