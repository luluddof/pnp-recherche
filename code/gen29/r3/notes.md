# G29-R3 adversarial review notes (start 14:16 UTC)
Targets: A = gen29/h1 (Majorana family, exact N=8, literature), B = gen29/x2 (composite SWAP->P_3, alpha*(P_3)).
## Plan
1 Majorana fact (own code, k=1..7) 2 family p,q (WHT k<=7; explicit operators N=8,16,32 mod 2 primes)
3 best-known balanced (HR, singlet, Hrubes finite-size?) 4 asymptotics 5 exact N=8 q<=4 own brute force 6 Servedio-Viola
7 X2 identity on exact examples 8 alpha*(P_3)=1/2
## Log
- 14:20 [E] maj_check.py (own code, RANDOM frames by backtracking, not JW): Majorana fact holds k=1..6 (all seeds) with AFFINE forms
  y_s = Q(s)+B(s,w); with pure LINEAR forms B(s,w) it FAILS (frames have mixed Q(s)) -> brief wording "linear forms" needs the
  Q(s) shift (H1's notes/proof do include it via c). y injective onto parity class: yes.
- fam.out: WHT counts k=2..6 reproduce H1 frontiers exactly with a different frame. q <= sum_{j<=d} C(2k+1,j) always
  (equality except 1 window); p = formula for every d>=1 window; for d=0 (constant f) p < formula (formula only an upper bound).
- expl8.out, expl16.out: explicit K = sum c(s) P_s(x)P_s, rank(SWAP-K) and realignment rank mod 1000003 & 998244353 equal
  (q, p from phi) for ALL windows at N=8,16 (incl. (8,8)@8, (10,46),(46,10)@16).
- 14:21 [E] bal_table.py/.out: balanced comparison incl. the project's Johnson-regraded Hrubes (JH, gen24/g24a, count n*min(C(t,j),N),
  t = rho(N)-1, all admissible k', j1+j2 = k'-1). k=18 (N=2^18): JH (j1=j2=2, k'=5, t=35) gives 2^18*595 = 155,975,680
  < Majorana 176,142,312.  => claim "best known for ALL 5<=k<=31 except k=4" REFUTED at k=18 (H1 never compared with JH).
  Majorana wins all other k in 5..31 vs JH/HR/N^2/8/singlet; also JH wins k=36; asymptotically JH 1.6147 < Majorana 1.6226.
- N=32 explicit (expl32.out): (12,232),(67,67),(232,12) CONFIRMED mod both primes.
- [E] fam7.out: k=7 (N=128) WHT with random frame: frontier identical to H1 incl. (576,576). 
- [E] exact8_mine.py (own: affine-orthogonal orbits of supports by union-find, vertex enumeration of G on image cells,
  sympy exact recheck): N=8 p_min(q)=28,28,28,24 for q=1..4 (4 orbits at q=3, 15 at q=4; optimum q=4 = one qubit {I,X,Z,Y});
  exact8_q5.py: q=5 -> 24 (110 O(Q)-orbits of {0}uT, float only) = H1. q=6 (=19, MILP) NOT checked by me.
- Lemma q<=2k => p>=N^2/4: proof re-derived (translate 0 into S, kernel of L nonzero, pairs {w,w+v} half split): valid.
  lemma_k2.out exhaustive N=4: p_min = 6,6,6,4,... for q=1..4 -> tight at q=2k (p=4=N^2/4).
- [N] asym.py: 2H(1/4)=1.62256; exact family balanced exponent k=100/400/1600/6400: 1.593/1.613/1.6195/1.6216 -> 1.6226 OK.
  Corner: delta*ln^2(1/eta)/eta^2 from the entropy formula = 0.073 (eta=.16), 0.16 (.023), 0.28 (3.6e-4), 0.42 (7e-9), 0.53 (1e-18)
  -> ln2 only as eta->0 (correction factor ln^2(1/eta)/ln^2(e/tau)). Leading asymptotics CONFIRMED; constant meaningless at real eta.
- Servedio-Viola (search snippets; eccc + columbia hosts blocked): highlighted special case = Hamming distance from truth tables of
  s-sparse polynomials over +-1 vars (= Bell model EXACTLY, s=q, fraction p/N^2); PROVED theorem = SIGN of s-monomial polynomial
  disagrees with IP on Omega(1/s): sign version (stronger than the exact-agreement line). Correspondence exact for the problem,
  not for the theorem; constants unspecified (Omega vs 1/(4q)).
- [E] p619.py: (q,p)=(6,19)@8 explicit operator check (osr 6, rank 19 mod 2 primes) CONFIRMED.
## B (X2)
- [E] x2_identity.py/.out: identity sigma - K12K23 = R12 S23 + K12 R23 holds exactly (random integer K, N=2,3,4, q=1..3; K=I N=2..5
  over Z; Bell (19,6)@8 and Majorana (8,8)@8 mod 1000003). sigma = S12 S23 maps (a1,a2,a3)->(a3,a1,a2) = P_3. Ranks: 4,16,40,80
  for K=I (= X2's), 288 <= 304 for (19,6)@8, 128 = 2pN (tight) for (8,8)@8. Flattenings of K12K23: site1|23 = q, site3|12 = q,
  site2|13 <= q^2: consistent with ring bonds (q,q,1) (and the explicit MPO sum_{l,m} G_l (x) G'_l H_m (x) H'_m is immediate).
  Algebra trivial: S12 S23 = (R12+K12)S23 = R12 S23 + K12(R23+K23). CONFIRMED.
- [N] alpha_p3_mine.py/.out (mpmath, gen27 formula): composite frontier alpha(mu) < 1/2 for all mu, h - ||.||_2 = (5/6)eps^2 > 0 at
  mu = 1-eps, alpha(1-eps) = .447,.465,.479,.487,.491,.495 for eps = 1e-2,1e-3,1e-5,1e-8,1e-12,1e-20 -> sup 1/2 not attained;
  max over mu<=1/2 is 0.173 (mu=1/2). Corner A: 1-u = (3/16)v^2 checked by hand. CONFIRMED (for the composite/known frontier,
  i.e. a lower bound on alpha*(P_3), design-level mapping assumed).
## VERDICTS (14:30 UTC)
A1 Majorana fact: CONFIRMED WITH CORRECTION (forms must be affine y_s = Q(s)+B(s,w); pure linear forms fail for generic frames)
A2 q,p formulas: CONFIRMED (q <= formula, equality generic; p = formula for d>=1, only an upper bound for d=0) N=8..128
A3 best known balanced 5<=k<=31: REFUTED at k=18 (JH/Hrubes count 155,975,680 < 176,142,312); true for the other k; 576@128 OK
A4 2H(1/4) and corner constant: CONFIRMED (asymptotic; corner constant converges extremely slowly)
A5 N=8: q=1..4 (and q=5) CONFIRMED independently; q=6 =19: upper bound confirmed, lower bound NOT CHECKED; lemma CONFIRMED
A6 Servedio-Viola: CONFIRMED WITH CORRECTION (problem = exact Bell model; proved theorem is the sign version, snippets only)
B1 identity: CONFIRMED.  B2 alpha*(P_3 composite)=1/2 as sup not attained: CONFIRMED (lower bound from known constructions)
