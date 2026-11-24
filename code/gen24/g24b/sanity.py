import math,numpy as np,sys
sys.argv=['x','0','A1']
exec(open('rate2.py').read().split("if MODE==")[0])
for (n,E0,H,exact) in [(3,3,4,8382),(3,2,4,3156),(2,4,8,8470),(2,3,6,1640),(1,6,12,862),(1,4,8,138)]:
    u=E0/n; v=H/(2*n)
    b=n*flat(u,v)+math.log(H+1)
    print(n,E0,H,"log exact %.3f  <=  log bound %.3f"%(math.log(exact),b), "OK" if math.log(exact)<=b+1e-9 else "VIOLATION")
