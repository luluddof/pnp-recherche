# G28-C3 notes (start 10:11 UTC). Routes H (Hagedorn) and P (projector chains) to a linear corner.
Tags: [P] proof sketch, [E] exact computation, [N] numeric, [C] counted only, [H] heuristic.
## 0. Plan (10:30)
- P1: "site-pattern model" (K = function of the per-site Sym/Anti pattern; = local-U(n)-invariant K on (F^n)^{(x)k}):
  random restriction mu_b = mixture of uniform subcubes + Smolensky over R  => delta <= C_n eta^2 (generalises Thm S to
  NON-symmetric functions). Verify formula for p,q vs explicit ranks; numeric search for non-symmetric g.
- P2: G-invariant model for groups with edge-type (geometric) dimension profiles (iterated wreath products / ultrametric).
- H1: symmetric orbifold Sym^M(Heis) dims in C1's GTC counting; H2 exp. many species / matrix (large-N) Fock spaces.
## 1. STOPPED at ~10:25 UTC (coordinator pause). No computations run yet. State of the analysis:
### 1a. Theorem S' (site-pattern / local-U(n)-invariant model is quadratic) [P, unreviewed, not yet numerically checked]
Setting: N = n^k, K in (x)_i span{1, SWAP_i} (equivalently K commutes with (U_1(x)..(x)U_k)^{(x)2} for all U_i in U(n)).
Then K = g(T), T in {0,1}^k = pattern of antisymmetric sites, trace measure mu_b, b = (n-1)/(2n); SWAP = (-1)^{|T|}.
Writing g = sum_U c_U prod_{i in U} tau_i (tau_i = SWAP_i - 1/n = centred variable), osr(K) = sum_{U: c_U != 0} (n^2-1)^{|U|}
(per-site left-factor space M_n = span{1} + traceless). Hence deg g <= log_{n^2-1} q.  p = mu_b(g != (-1)^{|T|}) N^2.
Proof of the bound: mu_b = mixture: each coordinate independently fixed to 0 w.p. 1-2b, else uniform on {0,1}.
On each random subcube (m free coords, uniform measure) g restricted has degree <= d and the target is the parity of the
free coords; Smolensky over R (y_S = parity * y_{S^c} on the agreement set) gives agreement <= 2^{-m} sum_{i<=(m+d)/2} C(m,i).
So p/N^2 >= E_{m~Bin(k,2b)} P(Bin(m,1/2) > (m+d)/2) ~ exp(-d^2/(2m)).  With d = eta k ln n/ln(n^2-1):
   delta <= eta^2 * n ln n / (2 (n-1) (ln(n^2-1))^2) + o(1)     (n=2: delta <= 0.574 eta^2; n large: -> eta^2/(8 ln n)).
Covers ALL functions of the singlet pattern (no permutation symmetry needed): generalises Thm S of G27-T3 (which was S_k-
symmetric Bell) for this submodel.  Does NOT cover the full Bell model (there a character costs 1, not 3^{|U|}; restriction
argument then only gives the row-uncertainty line bound B3).  Also extends (same proof) to O(n)-local-invariant K
(site algebra span{1,SWAP,Phi}) -- to be checked.
### 1b. Reflection-product model (Route P abstraction) [P/H]
SWAP = prod_i (1-2E_i), E_i commuting projectors; K = degree-d truncation in E's: p = mass(M>d), M = sum E_i,
q <= sum_{|S|<=d} osr(E_S).  Parity constraint P(M odd) ~ 1/2 is compatible with edge distributions (a0~a1~1/2, then
geometric).  T3's spectral test (|E[beta phi]| <= ||phi^||_1/N) excludes super-geometric tails but allows geometric tails
with rho >= 1/3 (Taylor truncation).  In Bell, Fourier-sparse indicators have granularity => osr(E_i) >= 1/eps_i; uniform
small eps design forced L ~ 2k/log k and is killed by the spectral test (coset indicators have ||^||_1 = 1).
### 1c. Route H [C, from reasoning only]
Sym^M(Heis) untwisted VOA: multisets of <= M Heisenberg states; for M >~ w growth log d(w) ~ c^2 w/(4 log w) (near-Hagedorn),
Cardy (c_eff = M) beyond; inside a VOA locality identity C1's span-doubling => endpoint excess => jump (expected; NOT computed).
Exp. many species / matrix large-N Fock spaces = Heisenberg with g(n) ~ e^{lambda n} generators: theta_end = 4 (C1 sect. 4).
### 1d. G-invariant model idea (not started): K invariant under diagonal conj. by G (e.g. iterated wreath S_n wr ... wr S_n,
whose orbitals on X x X have geometric = edge-type masses); p and q decompose over irreps -> small multiplicity matrices;
search block-vanishing patterns at N = 8, 9, 16.  Partition-algebra (G = S_N) case: trivial (p ~ N^2/2 or q ~ N^2 - O(N)) [P].
## 2. NEXT STEPS (precise)
 1. Script s_prime_check.py: for n=2, k<=5 build K = g(T) explicitly (L = sum SWAP_i), compare rank(SWAP-K), osr(K) (mod 2
    primes) with pattern formulas; then for k<=12 LP/L1 search of non-symmetric degree-d g under mu_{1/4}; check bound.
 2. Compute exact bound curve delta_bound(eta) at finite k vs singlet family (T3 bell_poly.py).
 3. Sym^M(Heis) dims -> C1 GTC curve (reuse gen28/c1/appb_curve.py logic).
 4. G-invariant wreath model at N=8,9,16.
Processes: none running.
## 3. Resumed 12:18 UTC (coordinator renamed my Theorem S' -> **Theorem S''**)
### 3a. S'' formula check [E] (s2_check.py, s2_check5.py)
K = sum_U c_U prod_{i in U}(SWAP_i - 1/2) built explicitly on (F^2)^{(x)k} (x)2, c_U solved mod p to agree with parity on random
pattern sets.  Pattern formulas p = sum_{T: g(T) != (-1)^|T|} 3^{k-|T|},  q = sum_{U: c_U != 0} 3^{|U|} == explicit
rank(SWAP-K), rank(realignment K) mod 2147483629 and 1000000007 in ALL cases k=2,3,4 (17 random instances, e.g. k=4:
(103,67),(70,88),(165,63),(100,28),(112,21)).  k=5 run: see s2_check5.out.
### 3b. S'' numerics [E small k / N heuristic]
- s2_exh4.py (exhaustive over all agreement sets, k=3,4): k=4 d=1 max mu-agreement 0.7734 > best symmetric 0.7383 (non-symmetric
  sets DO beat Thm-S symmetric ones at small k); d=2 .9492, d=3 .9961 (= bound!). All <= S'' bound (.8555,.9648,.9961).
- s2_greedy.py (randomised greedy, exact mod p consistency; reproduces exhaustive k=4):
  k=6: d=1 .6687 (sym .6526, bound .8062); d=2 .8591 (sym .8306, bound .9270); d>=3 = symmetric.
- s2_lp.py (L1 LP / IRLS) never beat the symmetric start: useless heuristic.
- s2_curve.py/.out: finite-k bound delta_bound(eta) = 0.574 eta^2 + O(log k / k) (slow convergence: delta/eta^2 at eta~.4:
  1.01 (k=40), .84 (k=80), .73 (k=160) -> .574).
- s2_greedy8.out (k=8, 60 restarts): d=1 .5961 (sym .5784), d=3 .9039 (sym .8862), d=4 .9735 (sym .9727), else = sym; all < bound.
  Non-symmetric optimisers = "d+1 full weight layers + partial next layer"; consistent with S'' (gain is O(1/sqrt k) only).
### 3c. G-invariant (permutation group) model, wreath products [E-ish: float ranks, block-vanishing search]
Model: K = sum_O c_O 1_O over orbits of G on X^4; p = rank(SWAP-K), q = rank(realignment). Pieces = central idempotents of
End_G(V(x)V) refined by SWAP (+/-). Search: random joint greedy over constraints "K = SWAP on piece" / "R(K) = 0 on piece",
generic point, true ranks (SVD).  wreath.py, wreath2.py, wreath_N8_9.out, wreath_N16_44.out, wreath_N16_2222.out, wreath_mult.out.
 S2wrS2wrS2 (N=8, 103 orbits): (28,1),(27,5),(26,9),(23,10),(22,12),(21,13),(18,14),(14,18),...,(1,28)   [self-dual]
 S4wrS2 / S2wrS4 (N=8): (27,6)/(27,5),(24,9),(22,12),(19,15) / (26,11),(18,14)
 S3wrS3 (N=9): (35,5),(34,7),(33,10),(28,15),(23,18),(21,21)
 S4wrS4 (N=16): (119,6),(116,9),(114,17),(112,29),(104,36),...,(56,69)  : p+q ~ N(N-1)/2 all along
 S2^{wr 4} (N=16, 227 orbits): (118,11),(115,13),(114,20),(108,21),(103,26),(93,34),(86,51),...,(60,68),(32,96)
 => NO corner: near the corner the saving is O(1)-O(N) for q ~ N^{0.6-1}; whole frontier ~ p+q ~ N^2/2 (far above the
 known general points, e.g. (22,5)@8 vs (27,5)@8, b(9)<=14 vs 21).  Multiplicities (wreath_mult.out): the big pieces are
 multiplicity-free (m=1, d = 24..54 at N=16), so K is a scalar there and block-vanishing is essentially exhaustive on the
 R side (partial ranks possible only in small m>=2 pieces).  Mechanism [H]: crossing (realignment) acts as a recoupling
 ("6j") transform between pieces of dimension d_rho >> 1; uncertainty is then ADDITIVE in piece dimensions (p+q >~ N^2/2),
 versus multiplicative (pq >~ N^2) for the abelian Bell case with 1-dim pieces.  Geometric (edge) orbital masses do not
 help: the Sym/Anti split inside each big piece is only available through SWAP itself.  Big non-abelian symmetry is BAD.
### 3d. Route H: symmetric orbifold Sym^M(Heis) [C] (symorb_count.py/.out, symorb2.py/.out; C1's appb_curve counting with
 p(w) -> d_M(w) = # multisets of <= M non-vacuum Heisenberg states)
 h=64: M=1 (App B) theta_end 2.69, jump T* 1.70, bal 1.86 | M=2: 2.53/1.39/1.72 | M=4: 2.53/1.31/1.69 |
       M=16: 3.11/1.67/1.84 | M>=64 (Hagedorn, ~ e^{c w/log w}): 3.26/1.70/1.85.
 Small M only rescales N (asymptotically = rank-M Heisenberg, C1 3a: subleading excess -> jump / circle); the Hagedorn
 regime (M >~ h) makes the endpoint WORSE (theta_end > 3 at h=64, -> 4 as h grows). No corner. As C1 predicted.
- s2_check_fast.py/.out: explicit permutation-expansion check, also qudits: n=2 k=5 (N=32, D=1024): (28,316),(430,12) OK;
  n=3 k=2: (9,17),(27,9),(45,8) OK; n=3 k=3 (N=27): (81,201),(108,89),(81,153) OK -- pattern formulas == explicit ranks
  mod two primes in every solvable instance (None = random interpolation system inconsistent, not a mismatch).
### 3e. Scope of S'' (when does the restriction+Smolensky method work?) [P sketch]
General site-local commutative model: K in (x)_i A_i, A_i commutative, SWAP_i in A_i, labels = minimal idempotents of A_i
(each in Sym or Anti), osr multiplicative over sites.  Method: write the per-site trace measure as a mixture containing,
with probability pi > 0, a uniform choice between one Sym label a and one Anti label b.  If every element of A_i that
separates a from b has osr >= c > 1, then deg (restricted g) <= log_c q and Smolensky gives
   delta <= (ln 2 / ...)  i.e.  delta <= eta^2 * (log2 N/k)^{-1} ... = O(eta^2)  [constant: d = eta log N / log c, m ~ pi k,
   error >= exp(-d^2/(2m))].
 * U(n)-local invariance (A = span{1,SWAP}): c = n^2-1 (proved above).
 * U(1)-local invariance for qubits (A = span{1, ZZ, XX+YY}; labels Phi(2-dim), Psi+, Psi-): pair (Psi+,Psi-), separating
   component XX+YY has osr 2 -> delta <= eta^2/ln 2 (bits) ~ 1.44 eta^2.  [same proof]
 * Full Bell (A = span{1,XX,YY,ZZ}): FAILS: XX (osr 1) already separates Psi- from Phi+ / Psi+; every Sym/Anti pair is
   separated by a single Pauli product of cost 1. So the Bell model is exactly the case this method cannot reach; there the
   question is low-rank Hadamard rigidity (T3 B1) -- still open.
- s2_greedy10.out (k=10, 6 restarts, d<=4): d=1 .5425 (sym .5318), d=4 .9291 (sym .9240); all far below bound (.75,.98).
- wreath_N16_2222_long.out (2500 samples, seed 7): (119,8),(113,16),(106,25),(100,34),(86,46),(69,65),(50,84),(40,88),
  (30,105),(1,120): p+q in [116,135] everywhere; q = N = 16 saves only 7 of 120.  Ranks by float SVD on 256x256 [N].
## 4. VERDICT (12:42 UTC)
Route P:
 * Site-local commutative models: Theorem S'' [P + formula checks E]: U(n)-local-invariant K (any function of the antisym
   site pattern, no symmetry) has delta <= eta^2 n ln n/(2(n-1) ln^2(n^2-1)) + o(1) (n=2: 0.574 eta^2); U(1)-local: ~1.44 eta^2.
   gamma = 1/2 at best.  Non-symmetric pattern functions beat symmetric ones at small k (exhaustive k=4, greedy k<=10) but
   only by O(1/sqrt k) mass.  Method provably cannot reach the full Bell model (cost-1 separators).
 * G-invariant permutation-group models (wreath products = ultrametric, geometric orbital masses; S_N/partition algebra):
   no corner; frontier ~ p+q ~ N^2/2 at N=8,9,16 [N]; big multiplicity-free pieces make crossing uncertainty additive.
 * Weyl qudit / U(x)U^dag families: (x)-commutation fails for odd n (T3), and {W_a (x) W_a^dag} generates a Heisenberg group whose
   commutant does not contain SWAP [P]: no commuting chain there.
Route H: Sym^M(Heis) [C]: small M = rescaled App B (gamma 1/2 / jump), Hagedorn regime M >~ h worsens the endpoint
 (theta_end ~3.3 at h=64, jump T*~1.7); exp. many species / matrix Fock = theta_end -> 4. No corner.
Best gamma: 1/2 (App B); nothing new realised.  Open: full Bell model (= low-rank Hadamard rigidity) and non-permutation
 abelian-invariant models (e.g. translation-invariant K, which contain Bell) -- not reachable by my methods.
Processes: none running (checked 12:42).
