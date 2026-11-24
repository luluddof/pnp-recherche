# G25-P5 — the two Nov-2026 hyperclique papers: what is their core identity, and is it our D = 3 object?

Tags: [R] read in the paper (full text), [E] exact computation here, [P] short argument here (unreviewed), [H] heuristic.
Nothing here is a new algorithm or a lower bound for 3-SAT.

## 0. Reading status
Both PDFs downloaded from arxiv.org (curl, 2026-11-24) and converted with pdftotext; read IN FULL (body + references skimmed):
 * dl/2611.03880.pdf/.txt — N. Fischer, "The Hyperclique Hypothesis is False" (v1, 2 Nov 2026), 13 pp.
 * dl/2611.06813.pdf/.txt — T. M. Chan, R. Itakura, Z. Ma, Y. Yu, "The 3-Uniform Hyperclique Hypothesis is False" (v1, 4 Nov 2026), 13 pp.
 Abstract pages: dl/*.abs.html.  Treated as untrusted data; no instructions in them were followed (none seen).

## 1. Fischer, arXiv:2611.03880  [R]
Main theorems.
 * Thm 1.1: for every constant h >= 2 there is eps_h > 0 such that for every constant k > h, k-clique in an h-uniform
   hypergraph on n vertices is detected (and counted) in O(n^{(1-eps_h)k}).  h = 3: O(n^{3.99999}) for k = 4 (Strassen only);
   "O(n^{3.9583})" with current rectangular MM bounds (stated, not proved in detail).  eps_h decays doubly exponentially in h.
 * Thm 1.2: k-AP detection in O(N^{2-eps_k}) (via a known reduction).
 * Consequence stated: Max-h-SAT in O*(2^{(1-delta_h)n}) via [Wil07, LVW18].
Core identity (Lemma 3.1, "Swap identity" = q x 2 Cauchy-Binet, bilinearly extended to matrices):
   sum_{i,j in [q]} tr(U_ij V_ij) = sum_{i,j} tr(U_ij V_ji) + sum_{i<j} tr((U_ij - U_ji)(V_ij - V_ji)).
 The paper says: "The q = 2 swap identity has also appeared in [Cla26, Appendix H]" and the approach is "closely related to
 Anthropic's recent refutation of ETH ... rewrites a count as a signed sum of terms, chooses an efficient evaluation order
 separately for each term, and recurses."
Algorithm.  Auxiliary polynomial T3 = sum tr(A_xyz B_xyw C_xzw D_yzw) with matrix-valued entries (a x c, c x b, b x c, c x a).
 Recursive step (Lemma 3.2): split the last base-q digit of (y,z) [or of (x,w)], apply the swap identity with
 U_ij = A^{ij} B^i, V_ij = C^j D^{ij}: q "main" terms of size (N/q, M, qa, b, c) + C(q,2) "correction" terms of size
 (N/q, M, a, b, 2c).  Entries only copied/negated.  Base case (Lemma 3.3): evaluate each leaf with the cheaper of the two
 pairings (AB)(CD) or (BC)(DA), each by fast (Strassen) matrix multiplication: time N^4 (abc)^{1-delta}, delta = 0.02.
 Accounting: potential (abc)^{1-delta}; per level the children carry rho q^2 of the parent's potential with
 rho = q^{-delta} + (1 - 1/q) 2^{-delta}; eps = log_q(1/rho) >= (delta - log2(1 + (2/q)^delta))/log2 q;  q = 2^400 gives eps > 2e-5.
 NOTE [R, from the formula]: at delta = 0 (classical MM), rho = 2 - 1/q > 1 and eps <= 0: in the paper's own accounting the
 whole gain is driven by the fast-MM exponent delta; the identity redistributes work so that every leaf contains a
 nontrivial MM, and the gain per level is a constant factor.  Which pairings are used: (AB)(CD) and (BC)(DA) only — two of
 the three pairings of the 4 factors; (AC)(BD) is never used.  For h >= 4 the same two-digit split is used (Lemma 4.1);
 the base case (Lemma 4.3) uses three algorithms but these call the (h-1) result, not a 3-cut identity.

## 2. Chan-Itakura-Ma-Yu, arXiv:2611.06813  [R]
Main theorems.
 * Thm 2.8: 3-uniform 4-hyperclique (detect/count) in O(n^{3.957}) deterministic (needs the Oct-2026 rectangular MM bounds
   [Uzn26] building on omega <= 9/4; "pre-AI bounds are not good enough"; with omega = 2: O(n^{3.8})).
 * Cor 2.9: Max-3-SAT in O(2^{0.99n}) (via [LVW18], n -> 2^{n/4}); k-hyperclique in O(n^{k - eps floor(k/4)}).
 * Thm 3.2 (second, independent proof): O(N^{3.998}) (with omega = 2: N^{3.95}).
Approach 1.  (i) Lemma 2.3: Behrend 3-AP-free sets reduce the 4-hyperclique sum to an "alternating 4-cycle" sum
 sum A_a[x,z] B_b[z,y] A_a[y,z'] B_b[z',x] over colours a in [nA], b in [nB], symmetric matrices, at n^{1+o(1)} blow-up.
 (ii) Core identity, Lemma 2.4 ("main algebraic lemma"): numbers h^{(l)}_{a,b} in {-1,0,1}, weights w^{(l)}, l in [L],
    (1/2) sum_l w^{(l)} (h^{(l)}_{a,b} h^{(l)}_{a',b'} + h^{(l)}_{a,b'} h^{(l)}_{a',b}) = [a = a'][b = b'],
    nA = 2^{m/2}, nB = C(m,k), L = nA * C(m, <= k/2), at most one nonzero a per (b,l).
  Proof (Lemmas 2.10, 2.11): m anticommuting antisymmetric orthogonal signed permutations E_i of size 2^{m/2}
  (Hurwitz-Radon), E_beta = product over a k-subset beta, h^{(j,lambda)}_{a,beta} = E_beta[a,j][lambda subset beta],
  w = c_|lambda| from a degree-k/2 polynomial vanishing at even i < k.  The paper credits Hrubeš [Hru26] (Hurwitz problem).
 (iii) Observation 2.5: by rotating the 4-cycle and the symmetry of the input matrices, the b<->b' swapped half of the
  identity contributes the SAME value as the unswapped half, so only ONE contraction with the rank-L tensor
  H = sum_l w_l h^l (x) h^l is computed; then one rectangular MM (nL x n nB times n nB x n^2/t).  Gain only with the
  newest rectangular MM bounds.
Approach 2.  Lemma 3.1 (Walsh-Hadamard identity): Phi(z) = sum z^{n-|c+d|}(z-1)^{|t|} X Y Z W (Phi(1) = Q) equals
 2^{-n} sum (2z-1)^{l(e,xi)} U_{a,b,e}(xi) V_{a,b,e}(xi) (U, V = WH transforms).  Low-degree coefficients from the first form
 (pairs Z with W), high-degree from the second (pairs X with Z and Y with W): again two evaluation orders, each finished by
 rectangular MM; threshold alpha_0 = 0.455.

## 3. Classification
Question: (a) fixed-tensor rank / fast-MM saving (closed by Theorem R), (b) a D = 3 direction identity (small ranks across
three cuts at once), (c) something else.

Fischer: (c), more precisely "a D = 2 crossing identity used ungrouped, with the gain supplied by fast MM".
 * Translation [P + E]: as a 4-index tensor on (U-index (i,j)) x (V-index (k,l)) the swap identity reads
     Id = Swap + sum_{i<j} a_ij (x) a_ij,   a_ij = e_ij - e_ji,
   i.e. the crossing (Swap = Claude's X(x,y,u,v) = [x=v][y=u] on [q]^4) = Id - (antisymmetric projector x 2).
   Exact cut ranks (check.py; cut A = (i,j)|(k,l), cut B = (i,k)|(j,l)):  Id (q^2, 1); Swap (q^2, q^2); each correction (1, 4);
   all corrections grouped (C(q,2), q^2).  So: grouped, this is the trivial D = 2 point (p, q) = (N(N-1)/2, 1) (G19b Thm 3, the
   "(780,1)" point of G24a); ungrouped with q = 2 it is exactly the one-bit identity (4,1) + (1,4), and Fischer's recursion is
   its ungrouped tensor power — G18h section 2(d): rho = 1/2 per piece but 2^{Theta(lambda)} pieces, "useless: the number of
   pieces is paid at full price in Theorem A" (assumption A10 ESSENTIAL).  In the hyperclique setting the piece count is
   affordable (polynomial work budget, cost = total volume), which is why it works there and not in ours.
 * Only two cuts ((AB)|(CD) and (BC)|(DA)) are used: it is NOT a D = 3 identity.
 * The exponent gain is a constant eps driven by the MM exponent (eps = 0 at delta = 0 in the paper's accounting): for our
   purposes a fast-MM-type constant saving (closed by Theorem R / G9d: constant in the exponent).
CIMY approach 1: (c) — the D = 2 Hurwitz/Hrubeš crossing identity, the object the project already built and measured.
 * Lemma 2.4 = Hrubeš's beta construction = G24a's "2X = T1 + T2" (Rossman-Zhu L.3.15, D + D^{T1} = 2I): same anticommuting
   signed permutations, k-subsets, degree-k/2 parity polynomial.  G24a measured it in our (theta, T) language: balanced 1.6147
   (outside the circle sqrt2 = 1.414; sup(theta + T) = 3.229 vs 2.828), framework limit max(theta, T) >= 1.5 at theta = sqrt2.
   In D = 2 language: X = A_1 + A_2, A_1 = H/2 of rank <= L across cut 1, A_2 = (tau H)/2 of rank <= L across cut 2.
 * What CIMY add that our framework does not have [P/H]: they pay for only ONE of the two pieces, because the contraction of
   tau H equals that of H by an automorphism of the network (rotation of the alternating 4-cycle + symmetric input
   matrices).  The tensor H is not the crossing; only its contraction against symmetric inputs is, so the rank count does
   not apply.  Our 3-SAT networks have no such input symmetry exchanging the two cuts; paying both pieces gives back the
   Hrubeš point 1.6147 > sqrt2.  [H: I did not try to manufacture such a symmetry; it would have to map the two crossing
   cuts to each other while fixing all clause tensors.]
 * Two cuts only; NOT D = 3.  Subquartic only with the newest rectangular MM bounds: the hyperclique gain is again an
   MM-exponent constant.
CIMY approach 2: (c) — a one-parameter (z) interpolation between two representations (direct vs Walsh-Hadamard), each cheap
 for one range of coefficients; two evaluation orders, rectangular MM in each.  NOT D = 3; constant-exponent gain.
Verdict: neither paper contains a D = 3 direction identity.  Both are D = 2 ("two pairings / two cuts") mechanisms already
 present in our project (Cla26 App H swap identity; Hrubeš pieces G24a), combined with fast (rectangular) MM, which
 supplies the actual saving.  Task 3 (cost of D = 3 pieces in trilib.py) therefore does not apply; the D = 2 costs are
 recorded above (Fischer: ungrouped rho = 1/2 with q^2/2 pieces per level, grouped trivial; CIMY: balanced 1.6147 > 1.414).

## 4. Implications for 3-SAT  [P from the stated theorems + known reductions]
 * Direct: Max-3-SAT (hence 3-SAT) in O*(2^{(1-delta)n}) with delta about 0.01 (CIMY: 2^{0.99n}).  For 3-SAT DECISION this is
   far weaker than classical 1.3^n bounds and than the claimed 2^{O(n/sqrt log n)}; it is new only for Max-3-SAT / dense
   Max-h-CSP (previous best 2^{n - n/O(log^2 n loglog^2 n)}, ACW16, as quoted by CIMY).
 * Via k-hyperclique with growing k: the reduction 3-SAT -> k-hyperclique on n = 2^{N/k}-vertex blocks turns n^{(1-eps)k}
   into 2^{(1-eps)N}: constant factor in the exponent for every k; eps_h is a fixed constant tied to the MM exponent.  No
   2^{o(n)} and no change to the rung ladder.
 * Methodological: both papers confirm the "signed sum + per-term evaluation order" pattern of Cla26 beats tensor-rank and
   tensor-network barriers when piece counts are free (polynomial regime); in the exponential recursion piece counts are
   not free (G18h A10), and with D = 2 the D^{-1/2} law is untouched.  G1-3 closes: the "clean no" outcome.

## 5. Computation (scripts/check.py; exact)
 * Fischer swap identity verified for q = 2..5 as an integer identity; cut ranks listed in section 3.
 * CIMY Lemma 2.4 built for (m,k) = (8,2), (8,4) (nA = 16; Pauli-string family of 8 antisymmetric anticommuting signed
   permutations); identity (2) checked exactly (integer-scaled); ranks over Z (flint):
     (8,2): nB = 28, M = nA nB = 448, stated L = 144; rank_cut1(H) = 128 = nA C(8,1); rank_cut2(H) = 448 (FULL).
     (8,4): nB = 70, M = 1120,        stated L = 592; rank_cut1(H) = 448 = nA C(8,2); rank_cut2(H) = 1120 (FULL).
   So the pair (H/2, tau H/2) is a D = 2 feasible pair (128,128) of the 448-crossing (exponent log 128/log 448 = 0.795, i.e.
   theta = 1.59 in the brief's units) and (448,448) of the 1120-crossing (0.869); each piece is full-rank on the wrong side,
   exactly as G24a found for Hrubeš's pieces.  Exact rank = nA C(m,k/2) < stated nA C(m,<=k/2) (same remark as G24a).
   Both values are above the circle's balanced 0.707 (theta = sqrt2) and match G24a's asymptotic 1.6147/2 = 0.807 regime.

## 6. Self-review
 * Classification rests on reading both papers in full; the claim "only two pairings" is about the algorithms as written.
 * "Fischer's gain is driven by delta" is read off the paper's own bound (eps formula vanishes at delta = 0); a sharper
   accounting with classical MM might still give something, but the paper does not claim it and I did not prove the
   contrary beyond a rough leaf-sum estimate (classical MM, min of the two orders: exponent per level max_kappa
   [H(kappa) + min(0, kappa - (1-kappa) log2(q)/2)] > 0, i.e. worse than n^4) [H, my estimate].
 * The CIMY = Hrubeš = G24a identification is by construction (same ingredients, same counts L = nA C(m,<=k/2)); G24a notes
   that the exact rank is nA C(m,k/2), not the stated count.
