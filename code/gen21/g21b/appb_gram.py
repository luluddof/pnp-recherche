# Appendix B over GF(p): rank mod p of the pairing (a,b) -> a[-h] b  on  V_{<=h} x V_h  (monomial bases, integer
# entries computed exactly with g20b's implementation of (47)), against p(h) (= rank over Q, Lemma 10) and against
# the number of partitions mu of h with z_mu != 0 mod p (no part divisible by p, all multiplicities < p).
import sys, math, importlib.util
spec = importlib.util.spec_from_file_location("ab", "/home/user/pnp-recherche/code/gen20/g20b/appb_exact.py")
ab = importlib.util.module_from_spec(spec); spec.loader.exec_module(ab)
from verify_gf import rank_mod
def gram(h):
    rows = [mu for w in range(h+1) for mu in ab.partitions(w)]
    cols = list(ab.partitions(h))
    M = []
    for mu in rows:
        M.append([dict(ab.apply_mono(mu, -h, nu)).get((), 0) for nu in cols])
    return rows, cols, M
if __name__ == "__main__":
    H = int(sys.argv[1]); primes = [2, 3, 5, 7]
    print("h  p(h) | rank_p Gram (p=2,3,5,7) | #{mu: z_mu != 0 mod p} | #{mu: no part divisible by p}")
    for h in range(1, H+1):
        rows, cols, M = gram(h)
        rk = [rank_mod(M, p) for p in primes]
        zc = [sum(1 for mu in cols if ab.zmu(mu) % p) for p in primes]
        nd = [sum(1 for mu in cols if all(m % p for m in mu)) for p in primes]
        print(f"{h:>2} {len(cols):>4} | {rk} | {zc} | {nd}", flush=True)
