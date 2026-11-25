# Independent R1 check. Level map of Lemma U (design level): e' = max_{(u,v) in H} min(e_a u, e_b v), H = down-closed
# convex hull of normalized piece costs. Computed GEOMETRICALLY: the max is where the ray (u,v) = s(1/e_a, 1/e_b) leaves H,
# so e' = s*. No support-function duality used.
import numpy as np, sys
def cr(o,a,b): return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
def convex_hull(pts):
    P=sorted(set(map(tuple,np.round(np.asarray(pts,float),14))))
    lo=[];up=[]
    for p in P:
        while len(lo)>=2 and cr(lo[-2],lo[-1],p)<=0: lo.pop()
        lo.append(p)
    for p in P[::-1]:
        while len(up)>=2 and cr(up[-2],up[-1],p)<=0: up.pop()
        up.append(p)
    return np.array(lo[:-1]+up[:-1])
class Hull:
    def __init__(self,pts):
        pts=np.asarray(pts,float)
        umax=pts[:,0].max(); vmax=pts[:,1].max()
        H=convex_hull(np.r_[pts,[[umax,0],[0,vmax],[0,0]]])
        H=H[(H[:,0]>1e-15)|(H[:,1]>1e-15)]
        ang=np.arctan2(H[:,1],H[:,0]); o=np.argsort(ang); self.V=H[o]; self.ang=ang[o]
    def s_star(self,ea,eb):
        # ray direction d=(1/ea,1/eb); find edge containing angle
        ea=np.atleast_1d(np.asarray(ea,float)); eb=np.atleast_1d(np.asarray(eb,float))
        d0=1/ea; d1=1/eb; a=np.arctan2(d1,d0)
        k=np.clip(np.searchsorted(self.ang,a)-1,0,len(self.V)-2)
        P0=self.V[k]; P1=self.V[k+1]
        # solve s*d = P0 + t(P1-P0)
        E=P1-P0
        det=d0*(-E[:,1]) - d1*(-E[:,0])
        s=(P0[:,0]*(-E[:,1]) - P0[:,1]*(-E[:,0]))/det
        return s
    def h(self,mu): return np.max(mu*self.V[:,0]+(1-mu)*self.V[:,1])
def fam(name):
    if name=='circle': t=np.linspace(0,np.pi/2,20001); return np.c_[np.cos(t),np.sin(t)]
    if name=='line': return [[1,0],[0,1]]
    if name=='p2q':  return [[1,0],[2/3,2/3],[0,1]]          # 2theta+T=4 and mirror
    if name.startswith('circ+lin'):   # circle arcs with phi in [phi0, pi/2-phi0] + corners (1,0),(0,1): finite corner slope
        phi0=float(name[8:]); t=np.linspace(phi0,np.pi/2-phi0,4001); return np.r_[np.c_[np.cos(t),np.sin(t)],[[1,0],[0,1]]]
    if name.startswith('lp'):
        p=float(name[2:]); t=np.linspace(0,np.pi/2,20001); c,s=np.cos(t),np.sin(t); r=(c**p+s**p)**(-1/p); return np.c_[r*c,r*s]
    raise ValueError
def fam_pow(g, s0=0.1):
    # circle with power corner at (1,0): v = k s^g, s = 1-u in [0,s0], joined to circle at s0 (one corner only)
    uj=1-s0; vj=np.sqrt(1-uj**2); k=vj/s0**g
    s=np.r_[0,np.logspace(-30,np.log10(s0),20000)]
    t=np.linspace(np.arctan2(vj,uj),np.pi/2,4001)
    return np.r_[np.c_[1-s,k*s**g],np.c_[np.cos(t),np.sin(t)]]
_old=fam
def fam(name):
    if name.startswith('pow'): return fam_pow(float(name[3:]))
    return _old(name)
