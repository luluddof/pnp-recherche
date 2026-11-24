# G9b (track B): is DECIDING cheaper than COUNTING in the 2026 frameworks?  (working notes, unreviewed)

Status tags: [P] proved here, [N] numerical only, [E] sketch, [O] open.

## 1. Ring-agnosticity of Claude's recursion
- Lemma 1 (network -> clique, balancing), Def. 1-4, Prop. 3 expansion: only + and x of table entries. Valid over any commutative ring. [P, by inspection]
- Only division: by d^{|R|} at each level (d = 2^t for Lemma 2; d = h! for App. B pieces, G4a step 2). Needs d invertible in the ring.
  * GF(2): eight-square pieces give 0 = 0 (d = 2^t = 0). Dead.
  * GF(p), p odd: Lemma 2 pieces fine (rho = 11/12).
  * App. B: needs p not dividing h!, with h ~ lambda^2 growing. So NOT available in any FIXED characteristic; available over GF(p) with p > h (p = polylog m).
- Characteristic 0 is used nowhere else (value bounds 2^{(2 lambda)^c} are only for bit length).
- Definition over a field F ("F-pieces with parameter rho"): for every lambda, poly(lambda) many A_i:( {0,1}^lambda )^4 -> F with
  sum_i A_i = c [x=v][y=u], c != 0 in F, and log2 rank_F(row flattening) + log2 rank_F(col flattening) <= 4 rho lambda + O(log lambda),
  factorizations computable in 2^{O(lambda)}. Then Theorem A(b) holds verbatim over F: Z mod char(F) (or any F-weighted Z) in 2^{O(m/(log m)^{log2(1/rho)})}. [E: I did not re-run Appendix A line by line.]
- Floor rho >= 1/2 - o(1) over EVERY field: the crossing is a permutation matrix of rank N^2 in both flattenings over any field. [P]
- Monotonicity: integer pieces with sum c X and p not dividing c reduce mod p with ranks <= Q-ranks. So rho_{GF(p)} <= rho_Q whenever p does not divide c. [P]
- Trace lemma [P]: both flattenings (indices (x,y),(v,u) and (x,u),(v,y)) have the SAME diagonal entries T(x,y,y,x). Hence tr_row(A) = tr_col(A) for
  every piece, and sum_i tr(A_i) = c N^2.  Two rank-one pieces that are both nilpotent force N^2 = 0 in F: possible only if char F divides N.
  This is the precise "characteristic-2 phenomenon" (analogue of sum-of-squares collapse): the vector g = e_01 + e_10 is isotropic over GF(2).

## 2. Small cases
One-bit (Pluecker) identity, c = 1, valid over Z hence every field:  [x=v][y=u] = [x=u][y=v] - g(x,y) g(u,v), g = [x != y](-1)^x.
So P_F(2) = 2 for every field (lower bound P >= N). SAT (satP.py): GF(2), N=2: (1 row,1 col) SAT; N=4: no decomposition with 4 terms (P_GF2(4) >= 5).
Tensor-power pieces A_i = [T^i](R + T C)^{(x)t}:
- Q: row rank = sum_{z<=i} C(t,z) 3^z; asymptotic rho = (2+log2 3)/4 = 0.8962 (paper App. H, N1).
- GF(3): R diagonalisable (g.g = 2 != 0); rank = sum over z<=i with 3 not dividing C(t-z,i-z); same asymptotic 0.8962. [P]
- GF(2): row matrix R = g g^T with g.g = 0: nilpotent, Jordan type J2 + 2 J1; col matrix of C the same.  A_i = unit * e_{t-i}(N_1..N_t) [P], so
  rank_row(A_i) = sum_n C(t,n) 2^{t-n} r(n,t-i), r(n,d) = sum_s rank_2 W_{s,s+d}(n) (inclusion matrices; Wilson 1990 formula, checked by brute force n<=10).
  Checked against direct 4^t x 4^t ranks for t <= 5 (direct.py).  Asymptotic rho_GF2 = 0.8113 (numerically = H(1/4)) [P modulo Wilson's theorem; limit evaluated numerically].
  finite t (rho_Q, rho_GF2): t=4: .758/.690; t=16: .837/.757; t=64: .875/.791; t=256: .889/.804.
- Exhaustive GF(2), lambda = 1: all 2^16 two-piece identities (exh2.py), 184 spectral class pairs; Pluecker is the unique best at t = 4, 5.
  CAVEAT: ranking at t=4,5 only; only the A_i grouping; asymptotics computed for Pluecker (0.811), R=0668 (0.906), R=0041 (1.0).
Table (asymptotic rho of best known pieces):
  lambda=1 tensor power:  Q 0.896 | GF(3) 0.896 | GF(2) 0.811
  eight-square (c = 2^t): Q 11/12 (0.901 by ranks) | GF(3) <= same | GF(2) unavailable
  App. B (c = h!):        Q 1/sqrt2 | GF(p), p > h: 1/sqrt2 | fixed p (incl. 2): unavailable
  floor:                  1/2 everywhere
So characteristic 2 is strictly better in the lambda=1 tensor-power class, but 0.811 > 1/sqrt 2: no gain in the exponent. [O]: F_2-pieces with rho < 1/sqrt 2.
Heuristic [E]: a tensor-power scheme built on a two-piece base identity of rank ~N per piece tends to rho -> 3/4, so beating 0.707 needs a non-tensor-power (graded, App. B-like) family in char 2.

## 3. Decision-only freedoms
(i) Valiant-Vazirani XORs: a dense parity is a 1-bit chain through all K nodes: +(K-1) edge bits per XOR, up to n XORs => + nK bits against S = 3m.
    Cost exponent e_l * nK with e_l ~ (log K)^{-1/2} >> 1/K: destroys the saving. Unnecessary anyway:
    - Z != 0 iff Z mod p != 0 w.h.p. for a random prime p in [n^2, 2n^2] (Z <= 2^n has <= n prime factors);
    - in char 2: Z(w) = sum_{sat a} prod w_i^{a_i} over GF(2^k), k = log n + O(1), random w (Schwartz-Zippel, multilinear); weights are unary tables, free.
    Both are "counting in another ring": same recursion, same rho. Only gain: poly-size arithmetic.
(ii) Hashing equality tests to lambda' < lambda bits: this relaxes the CSP (each cluster gets its own copy of shared variables, tied by lambda' parities).
    Over N: Z' = Z + junk >= Z: false POSITIVES (wrong side). Over GF(2): junk can also cancel: two-sided. A fixed bad difference survives w.p. 2^{-lambda'},
    but the set of realisable differences on an edge can be all of {0,1}^lambda, so soundness needs lambda' >= lambda. Even if sound it only rescales S (constant in the exponent).
    The sound direction (random affine restriction of edge strings) loses solutions: only helps when Z is huge, where sampling is better.
(iii) Monotone / Boolean semiring: pieces with OR instead of +, or nonnegative pieces. Each piece is then supported inside supp(X), a permutation
    matrix in both flattenings, so Boolean (and nonnegative) rank of a piece = its support size s_i in BOTH directions, and sum s_i >= N^2.
    Hence some piece has c_row + c_col >= 4 lambda - 2 log|I|: rho_bool = 1. [P] The method needs cancellation; "decision" enters only through a field.

## 4. Verdict
No decision-over-counting gain found. Every decision freedom is either a change of ring (same pieces, same rho) or unsound/costly.
A gain would follow from F-pieces with rho_F < 1/sqrt 2 over some finite field; char 2 shows a real small-case separation (0.811 < 0.896) but not below 0.707.
Context: classically decision 3-SAT (PPSZ, ~1.307^n) beats #3-SAT (Kutzkov 2007, 1.6423^n, quoted from memory, not re-verified) because PPSZ finds ONE solution; no such mechanism exists in a tensor-contraction framework.

## 5. Addenda / adversarial self-review
- Fourier remark on (i): in char != 2, [xor = b] = (1 + (-1)^b prod (-1)^{x_i})/2 is a sum of two unary-weighted networks: no edge bits but a factor 2 per XOR.
  So a VV parity costs min(1, e(K-1)) bits in the exponent; with up to n parities this is still too much. Not available in char 2 (needs 1/2).
- Framing: exact counting over Z itself reduces by CRT to GF(p) for ~n/log n primes p in [n, n^2]. So "counting" = large characteristic, and the ONLY decision-only
  algebraic freedom is a FIXED small characteristic (pieces that do not lift: c = 0 in char 0 lifts, nilpotent pieces). There the best known rho is worse (0.811 for p=2, 0.896 for odd p).
- Unproved / weak points: (a) Theorem A over F not re-derived line by line; (b) 0.8113 relies on Wilson's rank formula (checked n <= 10) and a numeric maximisation;
  (c) "GF(2) beats Q at lambda=1" compares against N1's search (R in {-1,0,1}^16, kappa in {1,2}, A_i grouping), not against all rational pieces;
  (d) exhaustive GF(2) ranking only at t = 4, 5 and only the A_i grouping (other symmetric or non-symmetric groupings untested);
  (e) lambda = 2 over GF(2): only P_GF2(4) >= 5 (SAT, Glucose, no proof log) and the t=2 tensor square; no exhaustive search (2^256 space);
  (f) GF(3) not searched exhaustively (3^16); (g) the 3/4 heuristic for tensor-power schemes is not a theorem; (h) OpenAI interpolation framework not analysed
  beyond noting that Hermite interpolation needs char > degree (division by factorials), so it also lives in large characteristic.
