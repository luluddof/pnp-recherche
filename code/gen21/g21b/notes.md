# G21b — does positive characteristic open the next step?  (working notes, saved as I go)
Tags: [P] short proof written here (unreviewed, same model family), [E] exact computation (own code), [N] numerical,
[L] literature from memory (not re-read), [O] open.  Nothing here is a new algorithm or a new bound for 3-SAT.
Setting: X = [x=v][y=u]; (p,q) feasible over K iff X = R + C, rr(R) <= p ((x,y)|(u,v)), rc(C) <= q ((x,u)|(y,v)).

## 1. Legitimacy: decision / counting mod p via Claude's recursion over F_{p^k}
Reduction.  3-CNF F on n variables.  Pick k with p^k > 2n (k = O(log n)), w in F_{p^k}^n uniformly; put the unary factor
w_i^{a_i} into one table containing variable i (tables of Definition 3 may take values in any commutative ring: Lemma 1,
Prop. 3 and Lemma 8 use + and x only, G9b [P, re-read here for Prop. 3 via Appendix A's proof: "each Z_iota is the value of the
row network and of the column network" — ring operations]).  Z(w) = sum_{a sat} prod w_i^{a_i} is multilinear, nonzero iff F is
satisfiable; Schwartz-Zippel: Pr[Z(w) = 0 | F sat] <= n/p^k < 1/2.  One-sided error, randomised.  Counting mod p: w = 1.
(Sparsification of Theorem A(d): decision needs only an OR over the sparsified formulas; counting mod p needs App. C's
disjoint version.)  Exact counting over Z is NOT obtained (it would need CRT over ~n/log n primes, i.e. large characteristic).
Steps of the paper that see the characteristic, and repairs:
 (S1) Division by d^{|R|} (Appendix A: "computing d^{|R|} and dividing by it").  Lemma 2 pieces: d = 2^t, dead in char 2,
      fine for p odd.  App. B: d = h! with h ~ lambda^2 -> h >= p for every fixed p: BREAKS.
      REPAIR [P + E]: Lemma 9(c) is an identity between integer matrices, sum_s Q_s(a,b) = (a[-h] b) Id, hence holds mod p.
      Replace Lemma 10 by a dual basis MOD p of the pairing V_{<=h} x V_h -> F_p, (a,b) -> a[-h]b.  Its rank mod p is
      d_p(h) := #{mu |- h : p does not divide z_mu} = #{partitions of h with no part divisible by p and all multiplicities < p}
      (out_gram.txt: equality checked for h <= 9, p = 2,3,5,7) [E].  With a_x, b_v dual mod p, y-states t^{nu_y}, u-functionals
      f_u = a_u[-h](.):  sum_s A_s = X mod p with d = 1, N = d_p(h).  Verified exactly mod p: (h,p) = (4,3),(6,3),(8,3),(6,5)
      (appb_modp.py) [E].  So App. B survives in every fixed characteristic, with N = d_p(h) instead of p(h).
 (S2) Integrality / bit sizes (Def. 7(c), "below 2^{(2 lambda)^c}"): irrelevant mod p (elements of F_{p^k} have O(log n) bits).
      Only simplifies.
 (S3) Lemma 10 (integral dual basis, h!/z_mu integers): breaks for h >= p; replaced as in (S1).
 (S4) Binomials C(k,j) for negative k in (47), Pascal/Vandermonde: integer identities, reduce mod p. Fine.
 (S5) Field size for isolation: F_{p^k}, k = O(log n), arithmetic poly(log n).  Decompositions found over F_{p^j} are used
      over F_{p^{jk'}} (take k a multiple of j); ranks of a FIXED matrix do not change under extension.
      CORRECTION to the brief: "rank invariant under extension => only p matters" is right for the conclusion, wrong for the
      reason.  The SET of decompositions grows with the field (as tensor rank does), so feasibility over F_p and over
      F_{p^j} can differ; what the algorithm can use is the frontier over the algebraic closure F_p-bar (any finite
      extension, uniformly computable), and even the BORDER frontier over F_p-bar (interpolation in eps needs only
      |K| > degree: extensions, G18g 1(e)).  An UNSAT answer over F_2 below excludes F_2-decompositions only.
 (S6) Monotonicity: an integer decomposition sum_i A_i = d X with p not dividing d reduces mod p with ranks <= the Q-ranks,
      so frontier_{F_p} <= frontier_Q on every p-integral point (G9b).  g20f/snf.py: all four known exact Q-points
      ((5,3)@4, (6,6)@5, (7,4)@5, (7,7)@6) have entries in {0,+-1} and are identities over Z: they hold over every F_p.
 Conclusion of 1: the decision / mod-p version is SOUND over F_{p^k}; the relevant object is the (border) frontier of the
 crossing over F_p-bar, plus the graded-profile analysis of App. B mod p (section 3a).

## 3a. Appendix B over F_p  [E for the identity and small ranks, P for the asymptotics]
 * Gram rank = d_p(h) (above).  Generating function prod_{p !| k} (1 - x^{pk})/(1 - x^k); by Meinardus,
   log d_p(h) ~ (1 - 1/p) pi sqrt(2h/3)  (p = 2: partitions into distinct odd parts, exp(pi sqrt(h/6)): tiny: d_2(9) = 2).
 * Exact ranks mod p of all thresholds (appb_modp_h*_p*.json) satisfy  rr(R_{s0}) <= D_p(2h+s0-1), rc(C_{s0}) <= D_p(2h-s0)
   with D_p(w) = sum_{w' <= w} d_p(w') (no violation in 4 runs), and are tight at the ends (rr = 1,2,4 at s0 = -15..-13, h=8, p=3):
   the pieces factor through the p-restricted quotient ("restricted Fock space": the radical of the form is a submodule,
   t_{pk} are never differentiated since nabla_{pk} = pk d/dt_{pk} = 0, multiplicities >= p die with the m! in z_mu).
 * Curve: theta(s), T(s) = log D_p(2h +- s) / log d_p(h) -> sqrt(2 +- sigma) (1-1/p) c / ((1-1/p) c) : the SAME circle.
   The factor (1 - 1/p) cancels between N and the costs.  rho stays 1/sqrt 2.  (If the ranks did not drop to D_p, the curve
   would be worse by 1/(1-1/p); they do drop.)  Finite N: balanced threshold value 7 at N = 4 (h=6,p=3), 29 at N = 6 (h=8,p=3),
   75 at N = 9 (h=6,p=5): like Q (G20b: App B trivial at small N), and with h larger for the same N.
 * What replaces "sum_s A_s = h! X": sum_s A_s = Gram(a_x,b_v) (x) Id_{y,u}, normalised by a mod-p dual basis (d = 1).

## 2. Exact small cases over GF(2), GF(3)  [E: sat2.py (CryptoMiniSat, native XOR), satp.py (CNF + CaDiCaL, one-hot GF(P)),
##    every SAT answer re-verified by verify_gf.py: own Gaussian elimination mod P + entrywise re-expansion R + C = X, R = UV]
 Calibration (agrees with known results): GF(2): (1,1)@2 SAT; (2,1)@3 UNSAT; (2,2)@3 UNSAT (= g19e q2n3.c exhaustive);
 (3,1)@3 SAT; (4,4)@4 SAT; (5,3)@4 SAT.  GF(3): (1,1)@2, (3,1)@3, (4,4)@4 SAT; (2,1)@3 UNSAT; (2,2)@3: CNF timed out (100 s),
 but excluded by G19b Thm 4 (char != 2).  Verified decompositions: dec_gf2_N*.json, dec_gf3_N*.json.
 Lower bounds valid over EVERY field (hence over F_p-bar): Theorem C (g19e, Lean) and Theorem 4 for N >= 4 (g20f: any char):
 q = 2 => p >= N.  Theorem 3 (q = 1: N(N-1)/2) needs char != 2; over GF(2) checked exhaustively N <= 4 (g20f).
 So for every field: b(2) = 1, b(3) = 3 (N = 3, char 2: (2,2) UNSAT over GF(2) only; over F_2-bar not proved), b(4) >= 3.
 Table (b = least p with (p,p) feasible; exact feasibility):
   N | Q exact                      | GF(2)                    | GF(3)
   2 | 1                            | 1                        | 1
   3 | 3                            | 3 (UNSAT (2,2))          | 3 (Thm 4)
   4 | 4 ((3,3) border only)        | 3 or 4: (3,3) see below  | 3 or 4: (3,3) not decided
   5 | <= 6 ((5,5) unresolved)      | <= 6 (Q point reduces)   | <= 6
   6 | <= 7                         | <= 7                     | <= 7

## 3b. Frobenius / Lucas  [P]
 * Frobenius: over F_{p^k}, R = UV  =>  R^{o p} (entrywise p-th power) = U^{o p} V^{o p}: rank preserved; same for rc.
   X^{o p} = X, and (R + C)^{o p} = R^{o p} + C^{o p} (freshman's dream: the mixed Hadamard terms carry C(p,j) = 0).
   So in char p the mixed terms DO cancel - but only for p-th powers of ONE decomposition, and the result is the
   Frobenius twist of the same decomposition, with the same (p,q).  A symmetry, not a gain.
 * Symmetric tensor powers: map every leg of X^{(x)p} (= crossing of size N^p) through the S_p-coinvariant quotient
   (or contract with symmetric tensors).  Then each S_p-orbit of mixed words R..C..R contributes C(p,j)(...) = 0, so
   pi(X^{(x)p}) = pi(R^{(x)p}) + pi(C^{(x)p}) with rr <= rr(R)^p, rc <= rc(C)^p: the mixed terms vanish.  But pi(X^{(x)p})
   has coefficient m(A) m(B) (multinomials) on the multiset pair (A,B), which is 0 mod p unless A, B are constant multisets:
   pi(X^{(x)p}) = X_N again (the Frobenius twist).  Sign-twisted legs (Lambda^p) give coefficient p! = 0: dead.
   With different decompositions in the p factors the words are not permutations of each other and nothing cancels.
   => "Mixed terms R(x)C are low rank in neither sense" is NOT avoided by characteristic p in any of these constructions.
 * Lucas: binomial matrices mod p are Kronecker products over base-p digits; a crossing decomposition built from them is
   digit-wise, i.e. a Kronecker / tensor-power construction, whose best known curves are G9b's (rho = 0.811 over F_2 from the
   Pluecker identity, 0.896 for odd p), above 1/sqrt2.  In App. B mod p the binomials C(l-1,n-1), C(-k-1,n-1) factor by Lucas
   but the operators also carry multiplicities (nabla_k = k mult_k): App. B mod p is not a pure Kronecker product, and its curve is
   the circle (3a).  No non-Kronecker recursive char-p family found.  [O]

## 3c. Non-rigidity over finite fields  [L from memory, P for the translation]
 * Dvir-Edelman (ToC 2019): p fixed, f: F_p^n -> F_p, M(x,y) = f(x+y) (N = p^n): for every eps there is delta > 0 with
   M = L + S, rank_{F_p} L <= N^{1-delta}, S with <= N^eps nonzeros per row (polynomial method / Croot-Lev-Pach: a polynomial
   of degree d in the 2n coordinates of (x,y) has rank <= 2 #{monomials of degree <= d/2}).  Dvir-Liu: the same over C for
   group-algebra matrices of abelian groups (f(x+y), Fourier); Alman-Williams: Walsh-Hadamard (f(x.y)) over Q.
 * Translation to the crossing.  Row flattening of X is delta(z + z') with z = (x,y), z' = (-v,-u): already a permutation, i.e.
   sparse; "rank + sparse" is vacuous for it.  What is needed is "rank + low operator-Schmidt rank", and sparsity does not bound
   the (x,u)|(y,v) rank (a permutation-sparse S can have rc = N^2).  So the rigidity statement itself cannot be fed in.
   What CAN be fed in is its proof technique: write X = prod_i (1 - (x_i - v_i)^{p-1})(1 - (y_i - u_i)^{p-1}) (a polynomial over
   F_p, char p essential) and group monomials by degree: every polynomial piece of degree D has rr, rc <= 2 #{monomials of
   degree <= D/2}.  This is a coordinate-wise (Kronecker-type) grouping, i.e. a sub-case of G9b's tensor-power pieces over F_p
   (rho >= 0.811 for p = 2 in that class), above the circle.  [P, no new curve computed]
 * Fourier side.  In characteristic p there is NO nontrivial p-th root of unity (t^p - 1 = (t-1)^p): the group algebra of
   F_p^n is local, there are no characters, no Weyl-Heisenberg basis for Z_p^n, and in char 2 the Pauli Z = diag(1,-1) = I.
   So the Pauli-diagonal dictionary of G20c (rank = support, osr = Fourier support) is DEAD in char 2 and for groups Z_p^n in
   char p.  It survives for Z_l^n with l != p (roots of unity in F_{p^k}, l | p^k - 1, 1/l^n exists), with the same
   statement as over C.  Exhaustive check, Pauli-diagonal family N = 4 (pauli_fp.py): least p for q = 1..6 is
   6,6,6,4,4,1 over GF(3) and GF(5) (running minimum of out_pauli_n2_p*.txt), identical to Q (G19b diag.py: 6,6,6,4,4).
   The C part costs Fourier sparsity (not rank) there exactly as over Q; the only char-p difference is that the finite-field
   uncertainty principle can be weaker (Evra-Kowalski-Lubotzky [L]); at N = 4 it makes no difference.  For odd l the
   transposition is not diagonal in the WH basis (G20c b'), also over F_{p^k}.

## 4. Verdict
 (iii) for the asymptotic question, and "no gain found" for the small cases:
 * Legitimacy: decision and #SAT mod p via Claude's recursion over F_{p^k} (k = O(log n), random weights, Schwartz-Zippel) is
   sound; the only broken step in a fixed characteristic is the division by h! (Lemma 10's dual basis), and it is repaired by a
   mod-p dual basis of the integer identity of Lemma 9(c) at the price N = d_p(h) (restricted partitions).  The relevant object
   is the border frontier over F_p-bar.
 * App. B mod p: same circle (the (1-1/p) factor cancels; exact ranks obey the restricted counts).  rho stays 1/sqrt2.
 * Every lower bound we have (line, Theorem C) is characteristic-free; Theorem 4 (q = 2 => p >= N, N >= 4) too.  Small exact
   cases: no pair found feasible over GF(2)/GF(3) and infeasible over Q (see section 2 for what was decided).  The known Q
   points are integral identities and hold over every F_p, so F_p is never worse at N <= 6.
 * Char-p mechanisms (Frobenius, symmetric powers, Lucas, Dvir-Edelman polynomial method) either reproduce the same
   decomposition (Frobenius twist) or are coordinate-wise Kronecker constructions (rho >= 0.811 in G9b's class).  The
   Fourier/Pauli family is dead in char 2 and identical to Q (N = 4) for odd p.
 * So: no strictly better asymptotic family; no ρ improvement for the decision problem.

## 2 (status of the open entries)  [machine shared with other agents: load 9-11 on 4 cores, each job got ~25-45% of a core]
 * (3,3)@4 over GF(2): split into the 4 cases rank G_1 = 1,2,3,4 (G_1 in rank normal form, legitimate by the GL_N x GL_N
   symmetry of X); every case hit the 590 s wall limit with NO answer (~2-5 min CPU each).  UNDECIDED.  Without symmetry
   breaking: no answer in ~8 min CPU.  (Over Q it is border-only, G19b.)
 * (5,2)@4 over GF(2): no answer in ~5 min CPU.  UNDECIDED.  (4,3)@4, N = 5 entries over GF(2): not run (CPU).
 * GF(3): the one-hot CNF is much weaker ((4,4)@4 SAT took > 100 s); (3,3)@4 over GF(3) not attempted.
 * So at N = 4 the only open question that could show a char-p small-case gain is whether (3,3) [or (4,3), (5,2)] is EXACTLY
   feasible over F_2 / F_3 (or an extension).  Even a YES would be a single small-N point: (3,3) is already border-feasible over
   Q and the border frontier is what the algorithm uses, and single points do not tensor (G19b sect. 6).
 * Nothing found: no pair feasible over GF(p) and infeasible over Q, no N with b_{GF(p)}(N) < b_Q(N), among the decided
   entries (N = 2, 3; and the q = 1 row for N <= 4 over GF(2), GF(3), g19e/g20f).

## Files
 notes.md; sat2.py (GF(2) SAT, CryptoMiniSat XOR); satp.py (GF(P) CNF, CaDiCaL); verify_gf.py (independent mod-P verifier);
 dec_gf2_N*.json, dec_gf3_N*.json (verified decompositions); appb_gram.py + out_gram.txt (Gram rank mod p = d_p(h));
 appb_modp.py + appb_modp_h*_p*.json, out_appb_h*_p*.txt (App. B with a mod-p dual basis: identity + exact threshold ranks);
 pauli_fp.py + out_pauli_n2_p*.txt (Pauli-diagonal family over GF(3), GF(5), N = 4).
 All my background jobs were stopped before the report.
