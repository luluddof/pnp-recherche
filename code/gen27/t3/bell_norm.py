# norm test: Q = sum |kappa| (Paulis unitary) for the polynomial family; B2 lemma: p >= (N^2 - N Q)/(1+Q)
import numpy as np
from fractions import Fraction
from math import comb, log2
from bell_poly import interp
from bell_poly_big import best_window
for k in [5,7,9,11]:
    N=2**k; n=N*N
    for d in range(0,4):
        A=best_window(k,d); f=interp(A)
        vals=[f(m) for m in range(k+1)]
        c=[]; diff=vals[:]
        for j in range(k+1):
            c.append(diff[0]); diff=[diff[i+1]-diff[i] for i in range(len(diff)-1)]
        # kappa at a character of weight u (support size u) has |coef| = |F(u)| where F(u)=sum_j c_j 4^-j C(k-u,j-u)
        Q=0.0; p=sum(comb(k,m)*3**(k-m) for m in range(k+1) if vals[m]!=(-1)**m)
        for u in range(k+1):
            F=sum(c[j]*Fraction(1,4**j)*comb(k-u,j-u) for j in range(u,k+1))
            Q+=abs(float(F))*comb(k,u)*3**u
        lb=(n-N*Q)/(1+Q)
        print("k=%d d=%d p=%d Q=%.3g  N^delta=%.3g  norm-lemma bound=%.3g"%(k,d,p,Q,n/p,lb))
