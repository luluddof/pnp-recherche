import sys; sys.dont_write_bytecode=True
from appb_exact import partitions, apply_mono
from math import comb
from fractions import Fraction
h=int(sys.argv[1]); mu=tuple(int(a) for a in sys.argv[2].split(',')); nu=tuple(int(a) for a in sys.argv[3].split(','))
sig=tuple(int(a) for a in sys.argv[4].split(',')); sigp=tuple(int(a) for a in sys.argv[5].split(','))
H=h+sum(mu); WM=4*h+2*sum(mu)+4
def ap(m, r, poly):
    out={}
    for s,c in poly.items():
        for s2,c2 in apply_mono(m, r, s): out[s2]=out.get(s2,0)+c*c2
    return out
fp={}; fm={}
for r in range(-WM, h+1):
    fp[-r-h]=ap(mu,r,ap(nu,-r,{sigp:1})).get(sig,0)
for r in range(-h, WM):
    fm[-r-h]=ap(nu,-r,ap(mu,r,{sigp:1})).get(sig,0)
def F(f,e): return sum((-1)**j*comb(H,j)*f.get(e-j,0) for j in range(H+1))
print("e  F+  F-")
for e in range(-2*h-2, H+3): print(e, F(fp,e), F(fm,e))
