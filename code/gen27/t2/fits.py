# trend fits on the frontier table: balanced exponent log_N b(N) = alpha + c/ln N; min r^2; min pq/N^2.
import math, numpy as np, json
import analysis as A
rows=json.load(open('trend_table.json'))
bE={r['N']:r['b'] for r in rows}
bST={9:14,10:15,11:20,12:22}          # g24d ST SAT arrays (exact, verified there; files not in repo)
bNUM={9:12,10:13}                    # G20b/G24e numerical (feasible-looking) values, lower end
best={N:min(bE[N],bST.get(N,10**9)) for N in bE}
def fitb(d,Ns):
    x=np.array([1/math.log(N) for N in Ns]); y=np.array([math.log(d[N])/math.log(N) for N in Ns])
    c=np.polyfit(x,y,1); return c[1],c[0]
print("N  b_exact(this table)  b_best(incl. ST)  log_N b  b/N")
for N in range(4,13): print(N,bE[N],best[N],f"{math.log(best[N])/math.log(N):.3f}",f"{best[N]/N:.2f}")
for lab,Ns in [('N=4..12',range(4,13)),('N=6..12',range(6,13)),('N=4..12 w/o 8',[n for n in range(4,13) if n!=8]),('N=8..12',range(8,13)),('N=9..12',range(9,13))]:
    a,c=fitb(best,list(Ns)); print(f"fit log_N b = alpha + c/lnN on {lab}: alpha = {a:.3f}, c = {c:.3f}")
# plain log-log slope
for lab,Ns in [('4..12',range(4,13)),('6..12',range(6,13)),('8..12',range(8,13))]:
    Ns=list(Ns); s=np.polyfit([math.log(N) for N in Ns],[math.log(best[N]) for N in Ns],1)[0]; print(f"log-log slope of b on {lab}: {s:.3f}")
