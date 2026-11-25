# counting proxy: dim U_W <= P(W) = sum_{w<=W} p(w); trivial-antidiagonal feasibility needs P(W)+P(W') >= N^2 = p(h)^2
from functools import lru_cache
@lru_cache(None)
def pp(n, m=None):
    if m is None: m = n
    if n == 0: return 1
    return sum(pp(n-k, k) for k in range(1, min(n, m)+1))
def p(n): return pp(n) if n >= 0 else 0
def P(W): return sum(p(w) for w in range(W+1)) if W >= 0 else 0
if __name__ == '__main__':
  for h in range(3, 41):
      N2 = p(h)**2
      best = None
      for W in range(0, 8*h):
          Wp = next(x for x in range(-1, 8*h) if min(P(W), N2) + min(P(x), N2) >= N2)
          k = 4*h-1-W-Wp
          if W <= Wp + 1 and (best is None or True): pass
      # balanced: W = W'
      Wb = next(W for W in range(0, 8*h) if 2*P(W) >= N2)
      Wb2 = next(W for W in range(0, 8*h) if P(W) + P(W-1) >= N2)   # W + W' = 2W-1
      kb = max(4*h-1-2*Wb, 4*h-1-(2*Wb2-1))
      print(f"h={h:>2} N={p(h):>6} balanced trivial W=W'={Wb:>3} (or {Wb2},{Wb2-1}): k_triv_bal={kb:>4}   k_triv/h={kb/h:+.2f}")
