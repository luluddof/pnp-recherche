# Idea 8: App B threshold decompositions reduced mod a small prime p > h (h! invertible, identity survives).
# Ranks mod p <= ranks over Q; does the frontier move?
import sys, math, time
sys.dont_write_bytecode = True
import flint
from appb_exact import pieces
h = int(sys.argv[1]); primes = [int(a) for a in sys.argv[2:]]
t0 = time.time(); N, mus, I, T = pieces(h); d = math.factorial(h)
print(f"h={h} N={N} pieces built {time.time()-t0:.1f}s", flush=True)
def flat_r(T4): return [T4[x][y][u][v] for x in range(N) for y in range(N) for u in range(N) for v in range(N)]
def flat_c(T4): return [T4[x][y][u][v] for x in range(N) for u in range(N) for y in range(N) for v in range(N)]
R = [[[[0]*N for _ in range(N)] for _ in range(N)] for _ in range(N)]
rows = []
for s0 in range(-2*h, 2*h+2):
    C = [[[[d*(x == v)*(y == u) - R[x][y][u][v] for v in range(N)] for u in range(N)] for y in range(N)] for x in range(N)]
    fr, fc = flat_r(R), flat_c(C)
    out = [s0]
    for p in [0] + primes:
        if p == 0:
            a = flint.fmpz_mat(N*N, N*N, fr).rank(); b = flint.fmpz_mat(N*N, N*N, fc).rank()
        else:
            a = flint.nmod_mat(N*N, N*N, [z % p for z in fr], p).rank(); b = flint.nmod_mat(N*N, N*N, [z % p for z in fc], p).rank()
        out.append((a, b))
    rows.append(out); print(out, flush=True)
    if s0 in T:
        for x in range(N):
            for y in range(N):
                for u in range(N):
                    for v in range(N): R[x][y][u][v] += T[s0][x][y][u][v]
for k, p in enumerate([0] + primes):
    print("p =", p if p else "Q", " balanced =", min(max(r[k+1]) for r in rows), " min rr+rc (nontrivial) =", min(sum(r[k+1]) for r in rows))
