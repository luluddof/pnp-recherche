# Independent brute-force check for n=3: enumerate all B2 circuits up to size 4, all 16 ops
import sys
n=3; N=8; FULL=255
ins=[sum(((t>>i)&1)<<t for t in range(N)) for i in range(n)]
def op(c,a,b):  # c in 0..15 truth table of op on (a,b)
    r=0
    if c&1: r|=(~a&~b)
    if c&2: r|=(~a&b)
    if c&4: r|=(a&~b)
    if c&8: r|=(a&b)
    return r&FULL
best=[99]*256
for f in [0,FULL]+ins+[FULL^x for x in ins]: best[f]=0
def rec(nodes,depth,maxd):
    for i in range(len(nodes)):
        for j in range(i+1,len(nodes)):
            a,b=nodes[i],nodes[j]
            for c in range(16):
                h=op(c,a,b)
                if best[h]>depth+1: best[h]=depth+1
                if depth+1<maxd: rec(nodes+[h],depth+1,maxd)
rec(ins,0,4)
from collections import Counter
print(Counter(best))
import json; json.dump(best,open('c3_brute.json','w'))
