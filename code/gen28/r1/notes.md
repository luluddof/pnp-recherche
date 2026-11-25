# G28-R1 review of G27-T3 (start 09:44 UTC)
## Plan / running notes
- Read notes + scripts. Writing independent code: singlet family built in computational basis from L = sum_i SWAP_i,
  K = f((k-L)/2) exact integers, rank mod 2 primes, osr = rank of realignment mod 2 primes.
- Lead (09:55): HR blocks [rho(N),N,N] tiled give integer SOS [N,N,N*ceil(N/rho(N))] -> (s,s) (RZ L3.15, G22b Thm 1).
  N=128: rho=16 -> s=1024 < 1156. Check.
- 09:50 singlet_check.py (computational basis, L=sum SWAP_i, exact ints, rank + realignment rank mod 1e9+7 and 998244353):
  N=8,16,32 all windows; N=64 best windows. All T3 points reproduced, incl. N=32 (16,349),(106,106),(376,16) [verify_k5_mine.txt].
  CORRECTION: closed forms p = sum_{m notin W}, q = sum_{u<=d} are only UPPER bounds: d=0 (constant f agrees on all m of same
  parity: SWAP-I has rank N(N-1)/2), window [0,k-1] (all strings of weight == k mod 2 vanish: q=28@8, 120@16, 496@32),
  k=6 W=[1,2] (identity coefficient vanishes: q=18 not 19). T3's code (bell_poly F(u) test) computes them exactly; the
  claimed points are right. Self-dual points (window [0,(k-1)/2]) exact by independent Krawtchouk formula (kraw.py) k<=15.
- kraw.py: HR tiling s = N*ceil(N/rho(N)) beats the self-dual point at k=7 (1024 < 1156). Family below HR at k=5,9,11,13,15.
- 10:05 hr_check.py: Hurwitz-Radon tiling [rho(N),N,N] blocks -> integer SOS [N,N,N*ceil(N/rho(N))]; s-square formula ->
  crossing pair (s,s) via D = sum vec(F)vec(F)^T, R = D*SWAP/2, C = SWAP - R (RZ L3.15; G22b Thm 1 STU reading).
  EXACT (two primes): (32,32)@16, (128,128)@32, (384,384)@64 [hr_k6.txt]; N=128: HR conditions + SOS identity verified,
  rank R <= rank D <= 1024 and osr C = rank D structurally -> (1024,1024)@128 < 1156.  Also (384,384)@64 < N^2/8=512 and <
  family's 694.  => "1156 at N=128 beats best known" REFUTED. N=32: 106 < 128 (HR = N^2/8 there) survives vs project's
  constructions; literature value of integer/complex [32,32,s] not checked (one web search: no table found).
  kraw.py table: family self-dual point beats HR tiling for k=5,9,11,13,15 (k=9: 12826 vs 14848), loses at k=7.
  Even k: family balanced worse than HR for k=6,8,10,12 (better at k=14).
- bell_mine.py: claim 1 (Weyl expansion + random q=N-1 truncations full rank N=3,4,5,8) and Bell formula == explicit rank
  (80/80 random cases N=4,8). bell_n4_mine.py: exhaustive N=4 Bell frontier q1-3:6, q4-5:4, q6-15:1 confirmed (both sides).
- corner.py / corner_lim.py: delta = (2/(3 ln2)) tau^2 exactly in the limit (two-sided window around 1/4, Bin(k,1/4)),
  eta = H(tau)+tau log2 3 -> delta ~ 0.96 eta^2/log2^2(1/eta) (very slow convergence: 0.31 at eta=1.6e-3, 0.43 at 2e-5).
  Finite k (300..3000) larger delta than limit (k^{-1/2} prefactor), still << eta^2/4. Claim 4 CONFIRMED.
- Claim 5 by hand: dimension proof correct, but needs a integer: p >= max_a (aN - q a^2) >= N^2/(4q) - q/4 (exact N^2/(4q) only
  when 2q | N). p_Bell(q) >= R_H(q): correct (Phi = H kappa^T H, rank <= nnz kappa; even R_H(rank kappa)); de Wolf N^2/(4r) for
  r <= N/2. Theorem S sketch: steps check (orbit-size bound via multinomial >= C(k,k-max); Krawtchouk polynomials of degree <= u;
  4 symmetric characters times beta are the 4 nontrivial characters of parity classes; tensor-grid unisolvence; sum over 8
  classes gives 0=8). Gaps to write: cube must fit in the type simplex (u <= ~k/8), mass estimate for u = Theta(k) is a KL
  bound exp(-k D) with D = Theta((u/k)^2), not a CLT one; fine for the stated conclusion.
