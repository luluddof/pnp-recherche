# own check of alpha*(P_3 composite) with the gen27/B2 formula alpha(mu)=sup{a: h(mu)<=||(mu,1-mu)||_{1/(1-a)}}, alpha*=sup_mu.
# frontier: image of App-B circle theta^2+T^2=4 under (u,v)=((1+theta)/3, 2T/3), truncated to u,v<=1, plus (1,0),(0,1).
import mpmath as mp
mp.mp.dps=50
def h(mu):
    mu=mp.mpf(mu)
    # max over phi in [0, phimax] of mu(1+2cos)/3+(1-mu)4 sin/3, sin<=3/4
    smax=mp.mpf(3)/4; phimax=mp.asin(smax)
    phis=mp.atan2(4*(1-mu),2*mu)                 # unconstrained maximiser
    ph=min(phis,phimax)
    val=mu*(1+2*mp.cos(ph))/3+(1-mu)*4*mp.sin(ph)/3
    return max(val,mu,1-mu)
def alpha(mu):
    H=h(mu); f=lambda a: (mp.mpf(mu)**(1/(1-a))+(1-mp.mpf(mu))**(1/(1-a)))**(1-a)
    if f(0)<H: return mp.mpf(0)
    lo,hi=mp.mpf(0),1-mp.mpf(10)**-40
    for _ in range(200):
        m=(lo+hi)/2
        if f(m)>=H: lo=m
        else: hi=m
    return lo
grid=[mp.mpf(i)/200 for i in range(1,200)]
vals=[(float(m),float(alpha(m))) for m in grid]
print('max alpha on grid mu in (0,1):',max(vals,key=lambda x:x[1]))
print('mu<=1/2 max alpha:',max((v for v in vals if v[0]<=0.5),key=lambda x:x[1]))
for e in [2,3,5,8,12,20,30]:
    mu=1-mp.mpf(10)**-e
    # margin at a=1/2: h - l2 norm (positive => alpha<1/2)
    print(f'eps=1e-{e}: alpha={mp.nstr(alpha(mu),8)}  h-||.||_2={mp.nstr(h(mu)-mp.sqrt(mu**2+(1-mu)**2),5)}')
