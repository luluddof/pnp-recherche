# G18b — charging the low-rank error part of the interpolation route differently

Agent notes (AI, same model family as the authors; not an independent human review; nothing here is formalised).
P = OpenAI preprint "The Exponential Time Hypothesis Is False" (itself unverified); read: §2.3–2.9, §3.1–3.3, §4.1–4.5.
L = ln n, L2 = lnln n. Status tags: [PROVED] = complete elementary argument written here (paper proof);
[SKETCH] = argument given, relies on lemmas of P not re-checked line by line; [NUMERIC]; [HEURISTIC]; [OPEN].

## 0. Summary
* The evaluator needs only the REACHABLE measurement tuples. Their number is between 2^rank and (D+1)^rank
  (both tight); P's box (2D+1)^rows is loose in P's own construction by a factor ~ r log2 r in the exponent.
* Charging the common error of P more cheaply (reachable tuples, or residues mod p) changes ONLY THE CONSTANT of
  n L2^2/sqrt(L): the binding constraint is rank fraction <= eps, i.e. R^2 >~ 2Q ln(1/eps), whatever the base.
* Changing BOTH the error structure (per-point errors, G7a's piecewise blocks, block length 3^a over F_3) and the
  charge (measurements mod 3, charge 3^rank) gives exponent O(n L2/sqrt(L)) — one loglog instead of two — but
  only for RANDOMISED DECISION (and for #F mod 3) of bounded-occurrence 3-CNF. [SKETCH, tiny cases verified]
* Exact counting: no gain beyond constants (the base is 2^{Theta(Q)} over Z, or n via CRT) — Theta(n L2^2/sqrt(L)).
* Lower bound in the framework with charge B^{M rho}: Omega(n sqrt(L2/L) min(ln B, L2)) and, for ln B >= L2,
  Omega(n L2^{1.5} sqrt(ln(ln B/L2))/sqrt(L)). So B = O(1): Omega(n sqrt(L2/L)) (G17c's estimate CONFIRMED as a
  lower-bound statement), B = polylog: Omega(n L2^{1.5}/sqrt L), B = n: Omega(n L2^2/sqrt L). For B = 3 the best
  upper bound is the O(n L2/sqrt L) above: gap sqrt(L2), OPEN.
* Verdict: P's interpolation bound IS improved by a factor lnln n in the exponent, in a restricted sense (decision,
  bounded occurrence). The headline is NOT beaten: on bounded-density inputs P has O(n/sqrt L); on unbounded-density
  k-SAT our variant only goes through sparsification and gets n (L2)^{O(1)}/L^{1/(2k)}.

## 1. What the evaluator needs; how many measurement vectors are reachable
Setting (P Def. 2.1/3.1, Lemma 3.2). At a node t: groups S_g, integer matrices C_g (c_g rows, row l1-norm <= D),
and factors A_g(x) that depend on x_h (h != g) only through y_h = C_h x_h. Then
    H(t) = sum_{x in {0,1}^n} prod_g A_g(x) = sum_{y in Y} prod_g B_g(y),  B_g(y) = sum_{u in {0,1}^{S_g}: C_g u = y_g} A_g(u; y),
where the sum may run over any set Y containing the reachable set prod_g Y_g, Y_g := C_g {0,1}^{S_g}
(unreachable y_g give an empty inner sum). P takes the box prod_g {-D..D}^{c_g}: (2D+1)^{rows}.

PROPOSITION 1 [PROVED]. C in Z^{k x m} of rank c; row i has positive part P_i, negative part N_i, P_i + N_i <= D.
 (a) For every row basis I:  2^c <= |C{0,1}^m| <= prod_{i in I} (P_i + N_i + 1) <= (D+1)^c;  also <= 2^{#nonzero columns}.
 (b) For a prime p: |C{0,1}^m mod p| <= p^{rank_p(C)} <= p^c.
 (c) Tight: c disjoint blocks of D ones (m >= cD) give exactly (D+1)^c; [I_c 0] gives 2^c; P's row (1,2,4,..)
     has c = 1, D = 2^m - 1 and 2^m values. With m columns and entries in {0,+-1} the maximum is between
     (floor(m/c)+1)^c and min(2^m, (m+1)^c).
Proof. Row i takes values in [-N_i, P_i] (P_i+N_i+1 values, not 2D+1). The other rows are fixed rational
combinations of the rows in I, so C x is determined by C_I x. Lower bound: choose c independent columns J; C is
injective on vectors supported on J. (b): the image lies in the column space over F_p.  QED
So: "a rank bound alone does not control the count" (P is right), but the count IS controlled by
rank x log(row norm + 1), and by rank x log p after reduction mod p. Nothing better than base
(row norm + 1) holds in general: the base is polynomial in the support size, not a constant.

In P's construction the box is loose in three independent ways [PROVED]:
 (i)  rows -> rank;
 (ii) (4.15) "keeps all rows, including zero rows": a basis vector psi = (x)_s psi_s restricted to S_{t,g}
      = {w : a_w(t) = g} is nonzero only if psi_t involves group g, and psi_t in {e_aa, e_ab +- e_ba} involves at
      most two groups. So at most 2 h c^B of the r h c^B rows are nonzero (r = 2: exactly c^B/2 rank per group,
      exp_aux_out.txt);
 (iii) |supp psi| <= 2^Q, so the row norms after restriction are <= 2^{Q-1}, not D = r^{2Q}.
Hence #reachable tuples <= (2^{Q-1}+1)^{2 h c^B}  against P's (2 r^{2Q}+1)^{r h c^B}: the exponent is smaller by
a factor ~ r log2 r. [HEURISTIC] The true count is 2^{Theta(Q) rank}: the rows are orthogonal with supports
~2^{Q(1-1/r)}, so y is close to a lattice Gaussian with ~Q/2 bits per coordinate. Tiny evidence: exp_reach_out.txt
(dense sub-boxes: log2 #reachable = 12.04 at rank 6, i.e. base 4, strictly above 2^rank; mod 3 exactly 3^rank).

## 2. Evaluators
E1 (exact, deterministic, over Z; "reachable box") [PROVED]. For each g compute a row basis I_g of C_{t,g} of small
norms (greedy) and the matrix T_g with C_{t,g} = T_g (C_{t,g})_{I_g}. Enumerate b = (b_g)_g, b_g in
prod_{i in I_g} [-N_i, P_i]; put y_g = T_g b_g; B_g = sum over u in {0,1}^{S_g} with (C_{t,g})_{I_g} u = b_g.
The classes {u : C_I u = b_g} partition {0,1}^{S_g} and y_g = C_{t,g} u on each class, so P's proof of (2.6)
applies verbatim. Polynomial space. Cost per restricted sum:
    prod_g prod_{i in I_g} (||row_i||_1 + 1) * sum_g 2^{|S_g|} * poly  <=  (2^{Q-1}+1)^{2hc^B} * r 2^{2n/r} * poly.
(Hashing the reachable set is not needed and would cost space.)

E2 (mod p; decision or #F mod p) [SKETCH; every step checked against P §2.4–2.5, tiny end-to-end test passes].
p odd prime, K = F_{p^k} with p^k > max(Q, 3n). Nodes tau_1..tau_Q distinct nonzero in K. Replacement rows over K
with P_{f,v}(0) = delta_v, degree <= R, and integer matrices C_{t,g} such that P_{f,v}(tau_t)|_{S_{t,g}} lies in
the K-row space of C_{t,g} for g != sigma_f(t). Weights w in K^n uniform;
    H_w(T) = sum_x prod_i w_i^{x_i} prod_f Phi_f(<P_{f,v}(T), x>),   H_w(0) = sum_{x sat} w^x.
 (1) Separation: if C_{t,h} x_h = C_{t,h} x'_h mod p then x_h - x'_h is in the K-kernel, so every factor assigned
     to g is a function of (x_g, (C_{t,h} x_h mod p)_{h != g}); the weight is local to groups. Lemma 3.2 is stated
     for a commutative ring; its proof ("every assignment has exactly one measurement tuple") is unchanged.
     Enumerate b_g in F_p^{rank_p C_{t,g}} for a row basis mod p (polynomial space):  p^{sum_g rank_p C_{t,g}} tuples.
 (2) Taylor data: (2.7)–(2.9) use only ring operations in K[z]/(z^J); no division by integers.
 (3) Recovery: H(tau_t + z) mod z^J for all t determines H when QJ > deg H (divisibility by prod (T-tau_t)^J),
     over any field; the confluent Vandermonde system is solved over K.
 (4) Decision: sum_{x sat} w^x is a nonzero multilinear polynomial in w iff F is satisfiable; Schwartz–Zippel:
     error <= n/|K| < 1/3, one-sided. With w = 1: H(0) = #F mod p, deterministic.
Cost per node:  sum_{j<=q} 2^j C(n,j) * p^{sum_g rank_p C_{t,g}} * r 2^{2n/r} * poly(n),  q = floor(deg H / Q).
Limits: E2 does not give #F. CRT needs primes up to ~0.7 n (product > 2^n), i.e. base n, worse than E1's base
2^{Q-1}+1 (= about D^{1/(2 log2 r)}). Z/p^e (Galois rings) needs the measurements mod p^e: base p^e, no gain.

P's own construction over F_p [PROVED]: needs p odd (the 1/2 in W_0 and the eigenbasis of the swap; the basis psi
has Gram determinant a power of 2) and p >= R (the R integers of I must be distinct mod p; f_I(j mod p) = (-1)^j
then holds on I). So E2 applies to P's rows with p the least prime >= R: charge p^{2hc^B}, ln p ~ L2/2.

## 3. Per-point errors in P's auxiliary setting, over Q and over F_p
P uses one common C^B for all nodes; Def. 3.1 allows C_{t,g} to depend on t. G7a's piecewise family transfers.
PROPOSITION 3 [PROVED; verified exactly for r = 2, Q = 3,4,5, p = 3 in exp_aux_out.txt]. K a field of
characteristic != 2, B = ([r]x[r])^Q, Pi_s the swaps, psi the eigenbasis (entries 0,+-1), z(psi) in {+-1}^Q its
sign pattern, |z| = #{s : z_s = -1}. Let m >= 1 and partition {0..Q} into intervals ("blocks") of <= m consecutive
integers. Assume for each block beta there are a^beta_0..a^beta_{m-1} in K with
        sum_{i<m} a^beta_i C(w,i) = (-1)^w   for all w in beta.                                   (#)
For |z| in beta put y_s(T) = [z_s = -1](1 - T/tau_s), F_beta(y) = sum_{i<m} a^beta_i e_i(y) (elementary symmetric),
        P_z(T) = 1 + (-1)^{|z|} ( F_beta(y(T)) - F_beta(y(0)) ),     P(T) = sum_psi P_{z(psi)}(T) psi psi^T/<psi,psi>.
Then (a) P(0) = I; (b) deg_T P <= m-1; (c) for every t: P(tau_t) = Pi_t A_t + E_t, where A_t acts as the identity
on position t (so row u of Pi_t A_t is supported on S_{t,d_u(t)}: "nice"), and the row space of E_t is spanned by
the psi^T with z(psi)_{-t} in Bad_t, Bad_t a subset of {z_{-t} : |z_{-t}| = top of a block, != Q}.
(d) (#) is solvable when char K = 0 or char K = p >= m (interpolation of degree < m at m points distinct mod p),
and when char K = p and m = p^a (Lucas: w -> C(w,i) mod p, i < p^a, is a basis of the p^a-periodic functions
Z -> F_p, and (-1)^w restricted to p^a consecutive integers extends to a p^a-periodic function).
(e) With W' ~ Bin(Q-1,(r-1)/2r), c_t := rank E_t <= D Pr[W' in block tops]; averaging over the m shifts of the
partition, some shift gives c_t <= D/m for all t simultaneously (the law of W' does not depend on t).
Proof. (a) F_beta(y(0)) = sum a_i C(|z|,i) = (-1)^{|z|} by (#), and P_z(0) = 1. (b) e_i(y(T)) has degree i.
(c) Pair z, z^t (z_t = +1, weights w, w+1): sum_{z_t} P_z Proj_z = Proj_{z_-t} [ (P_z+P_{z^t})/2 + Pi_t (P_z-P_{z^t})/2 ].
If w, w+1 lie in one block, then at T = tau_t both have the same y (coordinate t is 0 in both) and the same F, so
P_z + P_{z^t} = 2 + (-1)^w(F - (-1)^w) + (-1)^{w+1}(F - (-1)^{w+1}) = 0.  (e) is counting.  QED
Transfer to the variables: exactly P §4.5 with C_{t,g} := (diag(C_t^B,..,C_t^B) iota)|_{S_{t,g}}, zero rows dropped;
sum_g rank C_{t,g} <= 2 h c_t <= 2hD/m <= 4M/m <= 12 n/m (3-occurrence 3-CNF, M <= 3n). The support condition
(4.4) of P is not needed (it only served to bound the degree of P's phi_{u,y}). Lemma 4.5 (schedules, 12 D^2 <= M,
2Qr^2 <= n) is unchanged. End-to-end check of the separation with these C_{t,g} over F_9 on random tiny formulas
(n = 10, Q = 4, m = 3, weights, all nodes): exp_formula_out.txt, all OK.
Regimes: m << sqrt(Q): c_t/D ~ 1/m. m >> sqrt(Q) (one central block + tails): c_t/D ~ 2 sqrt(2/(pi Q)) e^{-m^2/2Q}
(P's common error: 2 e^{-m^2/2Q}; per-point errors gain sqrt(Q), G7a).

## 4. Re-optimisation and exponents (natural logs, divided by n; 3-occurrence 3-CNF, hD <= 2M <= 6n, 12 r^{4Q} <= n)
    E/n = 2 ln2 / r + 3 (R/Q) ln(2eQ/(3R)) + ERR,   deg P <= R,
 V0  P as written:              ERR = 6 r * 2e^{-R^2/8Q} * ln(2 r^{2Q}+1)
 V1  E1, common error:          ERR = 12 * 2e^{-R^2/2Q} * Q ln 2                     (Hoeffding tail)
 V1p E1, per-point:             ERR = 12 rho_pw(R) * Q ln 2
 V2  E2, p ~ R, common error:   ERR = 12 * 2e^{-R^2/2Q} * ln(2R)
 V3p E2, p ~ R, per-point:      ERR = 12 rho_pw(R) ln(2R)
 V3  E2, p = 3, per-point:      ERR = 12 ln 3 * rho_pw(R),  R = 3^a
Asymptotics. With a common error or with base 2^{Theta(Q)}, ERR <= eps forces R^2 >= 2Q ln(1/eps)(1+o(1)) with
ln(1/eps) = Theta(L2), while Q <= L/(4 ln r) and ln r >= ln(1/eps): R/Q = Theta(L2/sqrt L), Taylor term
Theta(L2^2/sqrt L). Leading constants 1.5 sqrt(2 kappa beta) (tail e^{-R^2/(kappa Q)}, R^2 = kappa Q beta L2):
V0 8.49 (kappa 8, beta 2), V1 ~3.7, V1p ~3.0, V2 2.12 (kappa 2, beta 1/2). ALL Theta(n L2^2/sqrt L). [PROVED as
a statement about these formulas]
V3: ERR = 12 ln3/R; minimise 12 ln3/R + 1.5 R ln Q/Q: R = sqrt(8 ln3 Q/ln Q), value 2 sqrt(18 ln3 ln Q/Q); with
r = sqrt(Q) (ln r = L2/2), Q = L/(2 L2):   E = (12.6 + o(1)) n L2/sqrt(L)  (x <= 1.16 for rounding R to 3^a).
V3p: ERR = 12 ln(2R)/R gives R ~ 2 sqrt Q and E = Theta(n L2^{1.5}/sqrt L): the Lucas trick (p = 3) is what
removes the last sqrt(L2).
[NUMERIC] opt_exponent_out.txt (continuous relaxation). E/n divided by L2/sqrt(L) for V3: 10.6, 11.3, 11.6, 11.8,
12.0, 12.2 at log2 n = 1e6, 1e9, 1e12, 1e15, 1e20, 1e30 (stable); optimal R/sqrt(Q) = 1.0 ... 0.37 (sub-sqrt(Q)
degree). V0, V1, V1p, V2 divided by L2^2/sqrt(L): slowly increasing towards their constants (e.g. V0 3.1 -> 6.5,
V2 1.47 -> 1.85). Ratio V0/V3 = 3.9, 7.3, 11, 15, 22, 36. All exponents are below n ln 2 only for log2 n >~ 1e5:
purely asymptotic statements.

THEOREM A [SKETCH: Prop. 1, 3 and E2 above are proved here; relies on P's Lemma 3.2, Lemma 3.3 / (2.9),
Lemma 4.5 and §4.5 being valid over F_{3^k}, which I checked from the overview proofs only; tiny cases verified].
Satisfiability of a 3-CNF with n variables, each in <= 3 clauses, is decided by a randomised one-sided-error
algorithm in time 2^{O(n lnln n/sqrt(ln n))} and polynomial space; #F mod 3 (any fixed odd prime) is computed
deterministically within the same bound. By the standard reductions the same holds for k-CNF with m <= Dn clauses
(constants depending on k, D).
Not obtained: exact counting; unbounded density. Via P §2.8 (sparsification + copies) decision of general k-SAT
costs O_k(sigma n + B n L2/sqrt L), B = (sigma^{-1} log(1/sigma))^{k-1}: n (L2)^{O(1)}/L^{1/(2k)}. P's better route for
unbounded density (§2.9: sparse F_2-systems, Beigel–Tarui modulo 2^{n+1}) is 2-adic and I found no odd-characteristic
substitute (the indicator of "L_i even" has degree ~ #monomials over F_3, multiplying deg H by polylog n). [OPEN]

## 5. Theorem 4 of G14d with the charge B^{M rho}  (verifying G17c's estimate)
Functional Ex_B := max( n/q, max_k rho_k M ln B, max_k hhat(h_k) ), hypotheses (H1), (H3), (H4) of G14d §4.
THEOREM 4_B [PROVED from G14d Theorem 2 and inequality (*), both taken as given; o(1) not explicit].
 (i)  Ex_B >= c n sqrt(L2/L) min(ln B, L2) for an absolute c > 0 and all B >= 2.
 (ii) If ln B/L2 -> infinity: Ex_B >= (1/2 - o(1)) u_1 n L2^{1.5}/sqrt(L), u_1 defined by
      c_0 e^{-u_1^2/2} u_1^{-3} ln B = L2/2 (c_0 of G14d Cor. 3), u_1 = (1+o(1)) sqrt(2 ln(ln B/L2)). For B = n this is G14d's
      (1/sqrt2) n L2^2/sqrt(L).
Proof. T := Ex_B = eps n; we may assume eps <= L^{-1/4}. q >= 1/eps, so ln q >= L2/4, and l <= (L+Lam)/(2 ln q - ln 2),
Lam = ln(M/n). x = d/N, u = x sqrt(l). By (H4) and monotonicity T >= n x ln(1/x) (x < 1/e, else T >= n/e).
By Theorem 2 and (*) with w = sqrt(l) (l >= 25; smaller l give rho >= 2^{-26}):
        max_k rho_k >= (1/(4 sqrt2)) l^{-1/2} exp(-(u+2)^2/2 - (u+2)^4/(6l)).
Case u <= 1: rho >= c_1/sqrt(l), T >= c_1 ln B M/sqrt(l) >= c_1 ln B n sqrt((2 ln q - ln 2)/L) (M/sqrt(ln M) is
increasing) >= c' ln B n sqrt(L2/L).
Case u > 1: if Lam <= L then x >= 1/sqrt(l) >= sqrt(ln q/(2L)) and T >= n x ln(1/x) >= c'' n sqrt(L2/L) L2.
If Lam > L, Cor. 3 and rho <= eps e^{-Lam}/ln B force u^2 >= 2 Lam - ln l - O(1), hence x^2 = u^2/l = Omega(ln q) > 1: absurd (l <= 4 Lam/ln q).
(ii): same with threshold u_1 instead of 1 (Cor. 3 of G14d): u < u_1 gives rho M ln B >= c_0 n ln B l^{-1/2} u_1^{-2}
e^{-u_1^2/2}(1-o(1)) = n u_1 (L2/2)/sqrt(l); u >= u_1 gives n (u_1/sqrt l) ln(sqrt(l)/u_1) = n u_1 (L2/2 - o(L2))/sqrt(l);
and 1/sqrt(l) >= sqrt(L2/L)(1-o(1)) because ln q >= (1/2 - o(1)) L2 when eps = L^{-1/2+o(1)}.  QED
Consequences: B = 2, 3: Omega(n sqrt(lnln n/ln n)) — G17c's unverified estimate is CORRECT as a lower bound of the
functional. B = polylog(n) (E2 with p ~ R): Omega(n L2^{1.5}/sqrt L). B = 2^{Theta(l)} (E1) or n: Omega(n L2^2/sqrt L).
Upper bounds in the same functional: B = 3: O(n L2/sqrt L) (piecewise blocks, §3–4) — gap sqrt(L2);
B = polylog: O(n L2^{1.5}/sqrt L) (V3p) — tight; B >= 2^{Theta(l)}: tight (G14d).
Where the gap for B = O(1) sits: at degree d = u sqrt(l) with u ~ 1/sqrt(ln l), Theorem 2 gives
rho >= ~0.4/sqrt(l) while the piecewise family has rho ~ 1/d = sqrt(ln l)/sqrt(l)... more precisely the optimum
uses d ~ sqrt(l/ln l) and rho ~ sqrt(ln l / l). Either a construction with rho(l,d) = O(1/sqrt l) at
d = O(sqrt(l)/ln l) exists (then 2^{O(n sqrt(L2/L))} for decision), or Theorem 1 of G14d is not tight for d < sqrt(l)
(it is certainly not tight at d = 0: truth 1, bound 1/sqrt(2 pi l)). [OPEN]
Remark: the sentence of G17c "the 2^{Omega(rank)} reading gives only Omega(n sqrt(loglog n/log n))" should carry the
qualifier that the best matching construction known gives O(n loglog n/sqrt(log n)), and only for decision / mod p.

## 6. Verdict
1. Does a different evaluator beat P's interpolation bound O(n L2^2/sqrt L)?
   - With P's rows (common error): NO. E1 and E2 divide the exponent by constants (asymptotically 8.49 -> 3.7 -> 2.12
     in the model of §4), the loglog power stays 2. [PROVED within the cost formulas]
   - With per-point piecewise rows + residues mod 3: YES for decision and for counting mod 3:
     O(n L2/sqrt L), loglog power 1. [SKETCH, Theorem A]. For exact counting: NO, power 2 (Theorem 4_B(ii) with
     ln B = Theta(Q) or L shows no accounting of the form B^{rank} can do better with such B).
2. Does it beat the headline? NO. Bounded density: P's App. D has O(n/sqrt L) (power 0) against our power 1.
   Unbounded density k-SAT: headline O(n L2/sqrt L) against our n (L2)^{O(1)}/L^{1/(2k)} (sparsification). Within the
   framework (charge B^{M rho}, B >= 2) the lower bound Omega(n sqrt(L2/L)) shows that even the best conceivable
   evaluator of this type cannot reach P's bounded-density bound n/sqrt L; reaching n sqrt(L2/L) is open.
3. Unproved / unverified steps:
   (U1) Theorem A relies on P's Lemmas 3.2, 3.3, 4.5 and §4.5; I checked that the overview proofs are
        characteristic-free but did not re-derive §3–4 in full, and P itself is not verified.
   (U2) Only r = 2, Q <= 5, m = 3, n = 10 were tested numerically (auxiliary kernel condition; evaluator identity).
        r >= 3 (equal pairs, sign probability (r-1)/2r) is covered by the proof only.
   (U3) The cost formulas of §4 are continuous relaxations with P's constants (hD <= 2M <= 6n); integrality
        (R = 3^a, prime p) costs constant factors; o(1) terms are not explicit; meaningful only for log2 n >~ 1e5.
   (U4) "#reachable = 2^{Theta(Q) rank}" for P's rows is heuristic; only 2^rank <= . <= (2^{Q-1}+1)^rank is proved.
        An exact-counting evaluator that does not enumerate measurement values is not excluded by anything here.
   (U5) Theorem 4_B inherits G14d's Theorem 2 and (*) (paper proofs; Lean covers only the weaker rank inequality).
   (U6) No odd-characteristic version of P's unbounded-density route was found; this is not a proof of impossibility.
   (U7) Gap sqrt(L2) between lower and upper bound for B = O(1).

## 7. Experiments (all exact, tiny; scripts in this directory)
* exp_aux.py / exp_aux_out.txt: Prop. 3 over F_3 (coefficients in F_9), r = 2, Q = 3,4,5, m = 3: P(0) = I, degree 2,
  Def. 3.1(4) kernel condition at every node with per-node C_t (rank over F_3), bad sets inside the predicted ones
  (some predicted-bad edges are accidentally good in characteristic 3), rank per group = c_t/2 (sum_g = c_t <= 2c_t).
* exp_formula.py / exp_formula_out.txt: 3 random 3-CNF (n = 10, 9 clauses, <= 3 occurrences, h = 1,2,3 copies),
  random weights in F_9: H_w(0) = weighted count; at every node the forms depend only on (x_g, y_{-g} mod 3) and
  sum over reachable y mod 3 of prod_g B_g(y) equals brute-force H_w(tau_t); #tuples <= 3^{sum rank_3}.
* exp_reach.py / exp_reach_out.txt: Prop. 1 tightness (blocks: exactly (D+1)^c; identity: 2^c); auxiliary rows on
  dense sub-boxes: 2^rank < #reachable <= basis-norm bound; mod 3: = 3^rank.
* opt_exponent.py / opt_exponent_out.txt: §4.

## 8. Adversarial self-review
* The gain comes from two changes together; neither alone helps: per-point errors with a large base (E1) stay at
  power 2 (G7a already said so), and a small base with the common error stays at power 2 (§4). Checked both.
* An earlier version of opt_exponent.py used min(1/R, central-block formula) for rho_pw; that is wrong for
  R < sqrt(Q) (the tails need their own blocks) and produced a spurious optimum at R = 1. Corrected to the mass of
  all block-boundary levels (~1/R); the numbers quoted are from the corrected script.
* Theorem A is a decision result with one-sided randomness; P's bounds are deterministic exact counting. Comparing
  loglog powers across these is not like for like; stated in the verdict.
* Possible objection: "the measurement classes mod 3 are coarser than the integer classes, is (2.5) still valid?"
  Yes: (2.5) is a linear identity over K and the kernel condition is over K; verified end to end (exp_formula).
* Possible objection: Hermite recovery with multiplicity J > p. It uses the truncated expansions H(tau+z) mod z^J
  (Hasse–Taylor), not derivatives; uniqueness is the divisibility argument, valid in any characteristic.
* Possible objection: in characteristic 3 the eigenvalue structure of W_0 collapses (j mod 3). True for P's
  f_I(W_0); Prop. 3 does not use W_0 but the symmetric functions e_i of the idempotents (1-Pi_s)/2, and Lucas'
  theorem gives period 3^a, which is exactly the block length. Verified for a = 1 only.
* The lower bound (ii) for B = n reproduces G14d's constant 1/sqrt2, a consistency check of §5.
* I did not check P's App. D, so "headline not beaten" compares with P's stated bounds only.
* Nothing here has been checked by a second party.
