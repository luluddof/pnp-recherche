import sys, fko
n=int(sys.argv[1]); a=float(sys.argv[2]); fko.CAP=10**7
m=int(a*n**1.5); r=fko.run(n,m,seed=1000*n+int(10*a),roots=(30 if n<=100 else 8),do_gauss=False)
print(' '.join(f'{k}={v:.4g}' if isinstance(v,float) else f'{k}={v}' for k,v in r.items()),flush=True)
