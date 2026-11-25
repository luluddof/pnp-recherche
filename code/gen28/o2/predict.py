# Check the spectral formula for merged ranks against merged_runs.out (exact, mod 2 primes).
#   row rank(A_i) = sum_{a+b+c=k, a<=i<=a+c} multinom * n0^a n1^b ng^c   (spectrum of M_row; i = #C factors)
#   col rank(A_i) = sum_{a'+b'+c'=k, a'<=k-i<=a'+c'} multinom * n0'^a' n1'^b' ng'^c'  (spectrum of M_col)
import re, sys, math
from pieces import *
from asym import spectrum_counts
def count(k, i, n0, n1, ng):
    t = 0
    for a in range(k+1):
        for b in range(k+1-a):
            c = k-a-b
            if a <= i <= a+c: t += math.factorial(k)//(math.factorial(a)*math.factorial(b)*math.factorial(c))*n0**a*n1**b*ng**c
    return t
txt = open('merged_runs.out').read().split('== ')[1:]
allok = True
for blk in txt:
    name = blk.split(':')[0]
    if name == 'oct': N, (R, C, s) = octonion()
    elif name.startswith('F1_'): N, (R, C, s) = F1(int(name[3:]))
    else: N, (R, C, s) = load_json(name)
    S = s*swap(N)
    cr = spectrum_counts(rowflat(R) @ np.linalg.inv(rowflat(S).astype(float)))
    cc = spectrum_counts(colflat(C) @ np.linalg.inv(colflat(S).astype(float)))
    rows = re.findall(r"k=(\d) i=(\d): merged row rank \[(\d+), (\d+)\]", blk); cols = re.findall(r"merged col rank \[(\d+), (\d+)\]", blk)
    for (k, i, r1, r2), (c1, c2) in zip(rows, cols):
        k, i = int(k), int(i)
        pr = count(k, i, *cr[:3]); pc = count(k, k-i, *cc[:3])
        ok = (pr == int(r1) == int(r2)) and (pc == int(c1) == int(c2)); allok &= ok
        if not ok: print("MISMATCH", name, k, i, r1, pr, c1, pc)
    print(f"{name.split('/')[-1]:22s} row spec (n0,n1,ng)={cr[:3]} vals {cr[5]} ss={cr[4]} | col {cc[:3]} vals {cc[5]} ss={cc[4]} | {len(rows)} classes checked")
print("ALL MATCH" if allok else "SOME MISMATCH")
