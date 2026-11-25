from symorb_count import symdims, curve
import math
for h in [40, 64]:
    for M in [1, 2, 4, h//4, h, 4*h+2]:
        d = symdims(4*h+2, M)
        N, te, jT, bal, pts = curve(d, h)
        print(f"h={h:3d} M={M:4d} log2N={math.log2(N):6.1f}  theta_end={te:.3f}  jump T*={jT:.3f}  balanced={bal:.3f}", flush=True)
