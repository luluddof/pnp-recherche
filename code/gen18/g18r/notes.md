# G18R — adversarial review of gen18/g18b (interpolation route of P, cheaper charge of the error part)

Reviewer: AI agent, same model family as the author of g18b; nothing here is formalised or human-checked.
P = OpenAI preprint "The Exponential Time Hypothesis Is False" (itself unverified; all verdicts below are
CONDITIONAL on P's framework §3–§4 being right, which I re-read but which nobody has verified independently).
Read in full for this review: P Def. 3.1, Lemma 3.2 + proof, §3.4 (Taylor/Hermite), Lemma 3.3 + proof, Prop. 3.4,
§4.1–4.6 (Lemmas 4.1–4.6, incl. 4.5 + proof and §4.5). All code here is mine (gf.py, design.py, ...), exact arithmetic.

## Verdicts
CLAIM 1 (Prop. 1: 2^c <= |C{0,1}^m| <= prod_{i in I}(||row_i||_1+1) for EVERY row basis I; mod p: <= p^{rank_p}) — CONFIRMED.
 * Lower bound holds for every integer matrix of rank c: take c independent columns J; C is injective on {0,1}^J.
   (C=(1,-1): 3 values, 2 <= 3 <= 3.) Upper: Cx is a function of C_I x and row i ranges over [-N_i,P_i]. Mod p: image lies in
   the F_p column space. The exact statement is "min over row bases"; (D+1)^c tight needs m >= cD.
 * claim1.py: 6000 random matrices (k<=4, m<=7, entries <=5, forced dependencies), 0 violations.
 * The consequences (ii) "<= 2 nonzero groups per row" and (iii) "norm <= 2^{Q-1}" for P's rows: re-derived, correct.

CLAIM 2 (E1, E2) — CONFIRMED (E1: elementary; E2: re-derived + tested end to end), with the precise hypotheses:
 * E2 does NOT run "P verbatim mod p": P's own rows cannot be reduced mod 3 (denominators (R-1)!Q!, nodes 1..Q collide
   mod 3, f_I(W_0) needs 1/j! and eigenvalues j as integers). E2 needs (a) nodes tau_t distinct, nonzero in K=F_{p^k}
   (|K| > Q), (b) rows over K[T] with P(0)=delta_v, (c) integer C_{t,g} with kernel condition over K. Given those:
   - Lemma 3.2 mod p: classes {u : C_g u = b mod p} partition {0,1}^{S_g}; differences lie in ker_K; proof unchanged. OK.
   - Lemma 3.3: Moebius inversion, signs only, any commutative ring. No factorial, no division, even for q >> p. OK.
   - Taylor data are the coefficients [z^j]H(tau_t+z) (computed in K[z]/z^J), i.e. Hasse derivatives, not d^j/dT^j. OK.
   - Recovery: P argues "rational Hermite matrix invertible over Q" (char 0). Replacement: CRT in K[T] for the pairwise
     coprime (T-tau_t)^J, deg H <= d_H < QJ. Valid in every characteristic. In e2e.py the confluent matrix
     [C(k,j) tau^{k-j}] was invertible over GF(9) in every run and the recovered degree was <= d_H.
   - Decision: sum_{x sat} w^x is multilinear, nonzero iff satisfiable; Schwartz–Zippel n/|K|; one-sided. OK.
   - Cost formula sum_{j<=q}2^j C(n,j) * p^{sum_g rank_p C_{t,g}} * r 2^{2n/r} * poly: correct.
 * Small inaccuracy: "P's construction over F_p needs p >= R" is right only with nodes moved to F_{p^k} (else p > Q).

CLAIM 3 (per-point blocks over F_3, E = (12.6+o(1)) n loglog n/sqrt(log n)) — CONFIRMED as a conditional theorem
 (conditional on P §3–4), with two corrections that do not change the bound:
 (C3-a) Prop. 3(d) is too restrictive and the "Lucas trick" is a red herring. det[C(lo+j,i)]_{j,i<m} = 1, so (#) has an
   INTEGER solution for every m and every block of m consecutive integers (Newton forward differences). The construction
   works over every field of characteristic != 2 for every m; m = 3^a is not needed (check_anym.py: m = 2,4,5,6,7 pass
   over GF(27)). The sentences "the Lucas trick (p=3) is what removes the last sqrt(L2)" and "x<=1.16 for rounding to
   3^a" are wrong attributions: what matters is (i) e_i of the idempotents instead of C(W_0,i) (no 1/i!), (ii) base 3.
   V3p (p ~ R) is a pointless variant: p = 3 is always available.
 (C3-b) the rank bound is loose by a factor >= 2: P_+ pi_+ + P_- pi_- = Pi_t(-P_-) + (P_+ + P_-) pi_+, so only the bad
   psi with psi_t symmetric (or only the antisymmetric ones: fraction (r-1)/2r) are needed as measurement rows
   (this is G7a's 1/(2(d+1))). Algebra only, not tested; it would lower 12.6 to about 8.9. Constant only.
 Re-derivation of the proof (own): P(T) = sum_z P_z(T) Proj_z with P_z(T) = (-1)^{|z|} F_beta(y(T)) (the notes' formula
 simplifies to this); P_z(0)=1 by (#); deg <= m-1; at tau_t, for z_t=+1 of weight w with w, w+1 in one block,
 P_z + P_{z^t} = 0, so the pair contributes Pi_t * (operator not touching position t): nice. Bad z_{-t}: w = top of a block;
 fraction Pr[Bin(Q-1,(r-1)/2r) in tops]; summing over the m shifts gives <= 1 in total, so some shift has <= 1/m for all t.
 Transfer (§4.5) and Lemma 4.5 are characteristic-free (see list). (4.4) is indeed not needed.
 Cost: E/n = 2ln2/r + 3((m-1)/Q) ln(2eQ/3(m-1)) + 12 ln3/m, constraint 12 r^{4Q} <= n. Own optimisation (opt.py, opt2.py,
 exact entropy form of ln C_q): E/n / (L2/sqrt L) = 10.2, 10.6, 11.1, 11.4, 11.6, 11.8, 12.0, 12.2, 12.3, 12.5 at
 log2 n = 1e5,1e6,1e8,1e10,1e12,1e15,1e20,1e30,1e50,1e100; optimum m = (2.9..3.0) sqrt(Q/ln Q), r ~ sqrt(Q);
 predicted limit sqrt2 * 2 sqrt(18 ln 3) = 12.58. Theta(n L2/sqrt L) inside this formula family is right (AM-GM on the
 last two terms with ln(Q/m) >= (1/2) ln Q, and ln r >= (1/2-o(1)) L2 forced by the 2ln2/r term).
 E/n < ln 2 only from log2 n ~ 1e5 on: purely asymptotic. P-style V0 in the same model: ratio V0/V3 = 1.8, 2.4, 3.5, 4.7
 at log2 n = 1e5, 1e6, 1e8, 1e10.
 Fairness of the comparison: NOT like for like, and the result is DOMINATED by P's own Theorem 1.1 (App. D, tensor
 contraction): deterministic exact #k-SAT in 2^{O(n L2/sqrt L)} for all k-CNF and 2^{O(n/sqrt L)} at bounded density.
 Theorem A (randomised decision / #F mod p, bounded density, O(n L2/sqrt L)) is weaker on every axis. What is improved
 is only the exponent attainable BY THE INTERPOLATION ROUTE (P's rho(n) = n L2^2/sqrt L, which P states for exact
 deterministic counting), and only for the weaker problems. It is a statement about the method.

CLAIM 4 (Theorem 4_B) — (i) CONFIRMED conditional on G14d Thm 2 and (*) (re-derived both cases u<=1, u>1, incl. Lam>L;
 the constant in "u^2 >= 2Lam - ln l - O(1)" is sloppy — the (u+2)^4/(6l) term costs a constant factor — conclusion
 unaffected). (ii) PARTLY: checked the two-case skeleton and the B = n consistency (1/sqrt2); o(1) terms and the
 Lam > L case not re-derived. Whether (H1)–(H4) of G14d (Rademacher model, p = 1/2) cover P's r-ary setting with
 p = (r-1)/2r is assumed by g18b, not shown. Consistency: Theorem 2 gives rho >= ~0.4/sqrt(l) at d << sqrt(l); the
 construction has rho ~ 1/m >= 1/sqrt(l): no contradiction; gap sqrt(L2) for B = O(1) is real and open.

## Characteristic-dependent steps of P (full proofs read) and what characteristic 3 does to them
 1. Def. 3.1(4): rows in Q[T], evaluation at integer nodes t in [Q].  Char 3: nodes collide -> nodes in F_{3^k}^x. FIXED by g18b.
 2. Def. 3.1(3) + Lemma 3.2: integer measurements in the box {-D..D}^c (integrality + archimedean bound).  Char 3:
    replaced by residues, 3^{rank_3}; proof of Lemma 3.2 unchanged. FINE.
 3. Def. 3.1(5), Lemma 4.6: denominators (R-1)!Q!, bit sizes (char 0 only). Char 3: Q! = 0, P's rows are not 3-integral:
    P's design itself does NOT reduce mod 3. Irrelevant once the rows are built over K directly.
 4. §3.4: "rational Hermite matrix invertible over Q; inverse acts by rational linear combinations" (A a Q-algebra,
    Prop. 3.4). Char 3: use CRT in K[T]; data are Hasse–Taylor coefficients already. FINE (tested).
 5. Lemma 3.3: none (signs only). FINE.  (3.14)–(3.15): counting only.
 6. (4.6) W_0 = (1/2) sum (I - Pi_s): 1/2. Char 3 fine, char 2 breaks.
 7. (4.7) f_I via binom(z-A_0, j) = (1/j!) prod: division by j! and integer eigenvalues j. Char 3: BREAKS for R > 3
    (j and j+3 have opposite parity but equal residue). g18b replaces it by e_i of the idempotents: FIXED.
 8. Lemma 4.1: orthogonal eigenbasis, <psi,psi> = 2^k; "orthogonal to retained rows => in span of the others".
    Char 3: norms are units, argument holds. Char 2: breaks (e_ab+e_ba = e_ab-e_ba).
 9. Lemma 4.2, 4.3, (4.12): linear algebra over any field.
10. (4.9) phi_{u,y} = prod (1 - T/s): division by s in [Q]. Char 3: breaks for 3 | s; not used by g18b (y_s = 1 - T/tau_s).
11. Lemma 4.5: combinatorial (labels), rational conditional expectations about labels only. Unaffected.
12. §4.5: iota, pi_c, C = C~ iota, "row sums <= D" (integrality, only for the box). Kernel verification is field-free.
13. §5–6 (Beigel–Tarui mod 2^{n+1}, sparsification): 2-adic / char 0; NOT used for bounded-occurrence clause products.
 Nothing in the path used by Theorem A breaks in characteristic 3.

## My numbers
* check_aux.py (matrix level, GF(9), all shifts, every node): r=2, Q=3,4,5 (rank test of Def. 3.1(4) with per-node
  C_t restricted to each group, non-vacuous kernels of dimension 14..110), r=3, Q=3 (same), r=3, Q=4 (D=6561, all 3 shifts) and r=2,
  Q=6 (D=4096, shifts 1 and 2 only: shift 3 hit my time limit) by support test of the nice part: P(0)=I, degree 2, nice part supported on a_y(t)=d_u(t), kernel
  condition holds; best-shift max_t c_t/D = 0.0625, 0.1406, 0.2109 (r=2, Q=3,4,5), 0.1111, 0.2222 (r=3, Q=3,4) <= 1/3.
  sum_g rank_3 C_{t,g} = c_t exactly (<= 2c_t). Some bad pairs are accidentally good in char 3.
* check_m9.py (eigen-pattern level, GF(27)): m=9 with Q=10,12,13 and m=3 with Q=12: all identities hold, degree 8 resp. 2;
  best-shift error fraction (r=2,3,4,8): m=9,Q=13: 0.0024, 0.011, 0.0114, 0.0147 (<= 1/9); m=3,Q=12: 0.327..0.333 (<= 1/3).
  m=9 at Q<=6 is the trivial single block. Matrix-level m=9 was not feasible (D >= 4^10).
* e2e.py: full evaluator (random labels -> copies -> rows over GF(9) -> per-node C_{t,g} mod 3 -> separated sums over
  reachable residue classes -> all (E,U) pairs of Lemma 3.3 -> confluent solve -> H(0)), checked against brute-force
  sum_{x sat} w^x (random weights) or #F mod 3, AND node by node against brute-force H(tau_t+z) mod z^J. See table below.
  The kernel condition was asserted on the actual rows for every (f,v,t,h != sigma); 'nonvac' counts cases where a
  class has >1 element and the row is non-constant on the group (i.e. the check could have failed).

### e2e.py runs (all PASSED: final value and every node's Taylor vector equal brute force; 0 failures)
 seed n cl r Q m shift | h | q  J | (E,U) pairs/node | rank_3 per group per node            | #classes        | nonvac kernel checks
  1   6  2 2 3 3  3    | 1 | 4  5 |   473            | [0,1],[1,0],[1,1]                    | <=2             | 3/18   (weighted)
  2   6  2 2 3 3  3    | 1 | 4  5 |   473            | [1,0],[1,1],[1,1]                    | <=2             | 3/18   (#F mod 3)
 11   8  3 2 3 3  3    | 2 | 6  7 |  5281            | [2,0],[1,1],[1,1]                    | <=4             | 6/27   (weighted)
 12   8  3 2 3 3  2    | 2 | 6  7 |  5281            | [3,1],[3,2],[4,1]                    | <=16            | 13/27  (#F mod 3)
 13   8  3 3 3 3  3    | 2 | 6  7 |  5281            | [0,3,1],[0,1,3],[0,2,0]              | <=8             | 4/54   (weighted, r=3)
 21  10  4 2 5 3  3    | 2 | 4  5 |  4521            | [4,1],[1,3],[2,4],[3,3],[4,1]        | <=16            | 18/60  (weighted, Q=5)
 14   8  3 2 3 3  1    | 9 | 6  7 |  5281            | [0,6],[6,0],[6,0] (one group empty)  | 64              | 18/27  (weighted, concentrated labels)
 15   9  4 2 3 3  3    | 9 | 4  5 |  2851            | all 0 (vacuous: rows of degree 1)    | 1               | 0/36   (#F mod 3)
 16   8  4 3 3 3  3    | 4 | 8  9 |  6561 (= 3^8)    | [0,0,4],[0,0,4],[0,2,4]              | <=16            | 0/72   (weighted, r=3)
 Not finished inside my time limit (killed, no result, NOT failures): n=10/Q=5/shift 1, n=12/Q=6, n=10/Q=4/5 clauses,
 n=8/Q=4. Limits of this test: n <= 10, m = 3 only, d_H small (2..5 clauses), so it tests identities, not the asymptotic
 regime; q < n was exercised (seeds 1,2,11,12,13,21,14,15) as well as q >= n (seed 16).

## Files (all in this directory)
 gf.py (GF(3^k) tables), design.py (block design), check_aux.py -> check_aux_small_out.txt, check_aux_big_out.txt;
 check_m9.py -> check_m9_out.txt; check_anym.py -> check_anym_out.txt; claim1.py -> claim1_out.txt;
 opt.py, opt2.py -> opt_out.txt (V0 column stops at 1e10: scipy binomial too slow beyond), opt2_out.txt (ignore the
 1e200 row: r-scan cap); e2e.py -> e2e_A_out.txt, e2e_B_out.txt, e2e_A2_out.txt (the tracebacks in e2e_A/B/C_out.txt are
 a bug of MY harness on an empty group, fixed before the *2 runs).

## Final judgement
There is a correct but narrow, conditional and non-headline result: IF P's §3–§4 framework is correct, then replacing P's
common error by per-point piecewise-block rows (any block length m, coefficients in F_{3^k}) and enumerating measurement
residues mod 3 gives a polynomial-space randomised one-sided decision procedure, and a deterministic #F mod 3 procedure,
for bounded-occurrence 3-CNF in time 2^{(12.6+o(1)) n lnln n/sqrt(ln n)} — one loglog factor below P's interpolation-route
exponent, but dominated by P's own Theorem 1.1 (App. D), so it improves the method's accounting, not the state of the problem.
