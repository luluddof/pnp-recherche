# G18e — closing the 4^u vs sqrt(u) gap in model DF (App. D tensor-contraction route)

Status: working notes. Level P = proof written here, not externally reviewed. Level E = evidence only.
(Interim save; sections are filled in as the work proceeds.)

## 1. Model DF and Theorem A in ten lines
1. Wires e = 1..M, each with two end bits (Y_e on a row vertex, X_e on a column vertex). I(Y,X) = prod_e [Y_e = X_e].
2. Directions i = 1..u: partitions pi_i of the vertices (App. D: group = symbol of the label in coordinate i, labels in [r]^a, r = 64).
3. A direction-i term is tau = prod_{g in pi_i} tau_g(end bits on the vertices of g), tau_g arbitrary.
4. DF scheme: I == sum_{k<=N} kappa_k tau^(k) entrywise mod 2^p, kappa_k in Frac(O); cost >= N.
   (App. D: k = (node t, zeta, mode tuple beta), direction v2(t); K_0 = I, (D.9), (D.16).)
5. Thm A: C1, C2 disjoint vertex sets, NO group of ANY pi_i meets both, w = #wires between C1 and C2 => N >= 2^w.
6. Proof: fix all other end bits to 0; I becomes a 2^w x 2^w permutation matrix P; each term becomes
   f(C1 bits) g(C2 bits), rank <= 1; P + 2^p E has det = +-1 mod 2, so rank 2^w over Frac(O).
7. Cor A1: C1 = prod S_i, C2 = prod ([r] \ S_i) with halves S_i: w = 2 * 4^-u * M.

WHERE IT LOSES.  The flattening is (C1 | C2) with C1, C2 required to be separated in every direction
simultaneously, so that every term has rank exactly <= 1.  That forces C1, C2 to be opposite ORTHANTS of
the label cube; each has density 2^-u, hence w = 2 * 2^-u * 2^-u * M.  The two factors 2^-u are:
(a) only packets whose symbol pair is split by S_i in all u coordinates survive (2^-u of the packets),
(b) inside such a packet only the two slots joining one antipodal pair of corners survive (2/2^u).
The rank used is the ordinary matrix rank over Frac(O), with the per-term bound "rank <= 1".
The loss is in insisting on rank <= 1 per term instead of "rank << 2^w per term" (see Thm A').

## 2. NEW: Theorem A' (unbalanced-group flattening), level P
Let C1, C2 be ANY disjoint vertex sets, w = #wires between them.  For a direction i and a group g of pi_i put
  e1(g) = #(C1-C2 wires whose C1-end lies in g),  e2(g) = #(C1-C2 wires whose C2-end lies in g).
Then   log2 N >= w - max_i sum_{g in pi_i} min(e1(g), e2(g)).
Proof.  Fix every end bit not on a C1-C2 wire to 0.  I restricts to a 2^w x 2^w permutation matrix P
(rows: bits at the C1 ends, columns: bits at the C2 ends).  A direction-i term restricts to the
Kronecker product over g of the matrices tau_g (rows: the e1(g) bits of g in C1, columns: the e2(g)
bits of g in C2), so its rank is prod_g rank(tau_g) <= 2^{sum_g min(e1(g), e2(g))} =: 2^{c_i}.
P + 2^p E = sum_k kappa_k tau^(k), and det(P + 2^p E) = +-1 mod 2O, so 2^w <= N * max_i 2^{c_i}.  QED
(Theorem A is the case c_i = 0.)

### Corollary A'1 (App. D's own labelling: 2^{Omega(M/sqrt u)}), level P
Take halves S_i of [r] in each of the u direction coordinates, s(label) = #{i : label_i in S_i},
C1 = {s > u/2}, C2 = {s < u/2} (u odd).  With uniform pair-class loads (Lemma D.3; the deviation of any
union of classes is <= M/sqrt(D_arr) = M/r^a by Cauchy-Schwarz, exactly as in Cor. A1):
 w = M/2,  e1(g) = w * mu_i^{C1}(g),  e2(g) = w * mu_i^{C2}(g)  (mu = i-th marginal of the uniform law on C1 / C2),
 sum_g min(e1,e2) = w (1 - TV(mu_i^{C1}, mu_i^{C2})) = w (1 - 2 eps),  eps = C(u-1,(u-1)/2)/2^u ~ 1/sqrt(2 pi u).
Hence  log2 N >= M * C(u-1,(u-1)/2)/2^u - O(M/r^a)  ~  0.399 M / sqrt(u).
So every DF identity on App. D's labelling has N >= 2^{(0.399-o(1)) M/sqrt u}, against App. D's
2^{O(M/sqrt u)} (if App. D is right): THE EXPONENTIAL GAP IN u IS CLOSED; in DF the truth is 2^{Theta(M/sqrt u)}.
Exact values (thmAprime.py, thmAprime_out.txt; fraction of M, uniform ordered-pair loads):
 u=1: 0.5 (= Thm A, tight), u=2: 0.25 (Thm A 0.125), u=3: 0.25 (Thm A 0.031), u=5: 0.1875, u=9: 0.137,
 u=41: 0.0627 (Thm A 4e-25), u=401: 0.0199.  Brute force over ALL cuts confirms 0.25 is the best this
 flattening family gives at u=2 and u=3.
Remarks.  (i) Unproved but expected: this family of flattenings cannot give more than O(M/sqrt u): the
 biases 2 eps_i are level-1 Fourier coefficients of the indicator of C1 (resp. C2), and
 sum_i eps_i^2 = O(1) by the level-1 inequality, so min_i eps_i = O(1/sqrt u).  Independently, App. D's
 upper bound (if right) caps every DF lower bound at O(M/sqrt u).
 (ii) General partitions (no product structure): with random signs sigma_i(g) and C1/C2 = {sum_i
 sigma_i(g_i(v)) >< 0} the same computation gives log2 N >= Omega(1/sqrt u) * min_i #(cut wires whose ends
 are separated by direction i); NOT written out in detail (level E).

## 3. Canonical instance and the Hurwitz link (level P for the reformulation, literature facts recalled, not re-derived)
One packet type, a = u: parties = corners of {0,1}^u, the wire of slot x joins corner x (row end) to the
antipodal corner ~x (column end); direction i cuts by bit i; EVERY wire is cut by EVERY direction.
With K = 2^{u-1} antipodal pairs of local dimension n_p = 2^{#wires on pair p}:
  N = min sum_i rank(B_i)  subject to  sum_i B_i^{Gamma_{S_i}} = Id  on  tensor_p C^{n_p},
Gamma_S = partial transpose on the pairs in S, S_i = {pairs whose orientation in direction i differs from direction 1}.
u = 2: Id = A + B^Gamma on C^n (x) C^n, equivalently "transpose on M_n = (map of Kraus rank R2) + (map of
rank R1)".  n <= N(n) (Thm A), and sigma_bil(n) <= N(n) <= 4 sigma_bil(n), sigma_bil = least k with
(sum x_i^2)(sum y_j^2) = sum_{l<=k} z_l z'_l, z, z' bilinear: the Hurwitz problem in the form of
Hrubes-Wigderson-Yehudayoff.  Known there: n <= sigma <= O(n^2/log n) (Hurwitz-Radon), and a bound
n^{1+eps} over C is a famous open problem (it would give exponential lower bounds for non-commutative
circuits).  So on a SINGLE packet type, improving Thm A's exponent M/2 (u = 2) is HWY-hard; the gain of
Thm A' comes entirely from mixing packet types (pair classes), where no such barrier applies.

## 4. Small cases (canonical single packet, r = 2 groups per direction)
Upper bounds: explicit integer decompositions, verified exactly (exact.py, exact_out.txt): terms summed
over Z equal I, and every term is checked to be rank 1 across the two groups of its direction.
Construction: Id = (Id - SWAP_{XY}) [direction A, n(n-1)/2 terms] + SWAP_{XY} [direction B, ONE term],
times the identity on spectator wires.  Lower bounds: Theorem A only (A' gives nothing more on a single
packet type).  "num" = real Levenberg-Marquardt search (search.py, out_*.txt): EVIDENCE, NOT PROOF.

 u  instance (slots)            M  trivial 2^M  proved LB  exact UB (verified)  numerics
 1  a=1, both slots             2      4           4            4               -            (exact min 4)
 2  a=2, slots {00,01}          2      4           2            2               -            (exact min 2)
 2  a=2, slots {00,11}          2      4           4            4               -            (exact min 4)
 2  a=2, slots {00,11,01}       3      8           4            4               N=3 not found (exact min 4)
 2  a=2, full packet            4     16           4            7               N=4, N=5 not found (30 starts each;
                                                                                 N=6 not run to completion)
 3  a=3, 3 slots on 3 pairs     3      8           2            4               N=2,3 not found; N=4 found
 3  a=3, 4 slots on 4 pairs     4     16           2            7               N<=4 not found; N=5 profile (3,1,1):
                                                                                 residual^2 -> 1.3e-10 with parameters
                                                                                 growing (164 -> 422): a BORDER-type
                                                                                 approximation, NOT exactly verified;
                                                                                 two-direction profile (3,2,0) fails.
Reading.  (i) App. D's count is far above 2^M at these sizes (rough estimate log2 #tuples ~ 44 M/sqrt(u),
 from p_t ~ 24 M D_t/u and 6M/D_t modes per weight; not verified), so no comparison with App. D is
 meaningful below u ~ 10^3.  (ii) log2(N)/M: u=1: 1; u=2: 0.5 (one wire per pair) and in [0.5, 0.70]
 (full packet); u=3: 0.667 (3 pairs, numerical) and ~0.58 (4 pairs, numerical border value 5; proved only
 [0.25, 0.70]).  These points do not follow 1/u, and are too few and too small to confirm 1/sqrt(u);
 they show (a) Thm A is NOT tight on single packets (numerically), (b) a third direction does help
 (5 < 7), apparently only through border-type (unbounded-coefficient) identities, which DF allows.
 (iii) The law itself is settled at the level of the model by Thm A'/Cor A'1, not by these cases:
 on uniformly loaded product labels the exponent is >= 0.399 M/sqrt(u), so NO DF scheme reaches M/u.
 Exact A' values on the uniform 2-colour instance: u=2: M/4, u=3: M/4 (trivial upper bound M/2).

## 5. Verdict
 - In DF on App. D's labelling the truth is 2^{Theta(M/sqrt u)}: lower bound proved here (Cor A'1,
   constant 0.399), upper bound = App. D taken at face value (constant of order tens).  The exponential
   gap in u is closed from the lower-bound side.  The sister agent's heuristic law "gain = (number of
   label directions)^{-1/2}" is, inside DF, a THEOREM as a lower bound: pieces that factor over the groups
   of ONE direction cannot reach the 1/u floor.
 - Since u <= a/2 and (D.6) forces r^{4a} <= M, u = O(log M): every one-level DF scheme with this packing
   costs 2^{Omega(M/sqrt(log M))}.  No scheme better than App. D was found; nothing changes for the
   3-SAT exponent (2^{O(n/sqrt(log n))} at best by this route, constants aside).
 - What could still beat sqrt(log): (1) many more directions than log M (A' only gives M/sqrt(u); whether
   DF identities with 2^{O(M/sqrt u)} terms exist for u >> log M is OPEN - App. D's packing (D.6) is what
   limits u); (2) terms that do not factor over one direction's groups (outside DF); (3) multi-level
   effects (A' counts terms at one level only).

## 6. Adversarial self-review
 - A' is a five-line flattening argument; checked: Kronecker structure of a restricted term (each cut-wire
   end lies in exactly one group), random exact rank test on a toy (thmAprime.py part 1), unimodularity
   mod 2 as in Thm A.  Not externally reviewed.
 - Cor A'1 uses Lemma D.3's uniform pair loads at face value, group = symbol in a direction coordinate
   for row AND column vertices, and halves of [r] (r even).  Error term (2r+1) M / r^a not optimised.
 - It is a bound on the NUMBER OF TERMS of a one-level identity, not on running time: children could be
   cheap, and schemes outside DF are untouched.
 - Upper bound side: nothing proved beyond tiny exact constructions.  App. D's constant is my rough estimate.
 - Small-case minima: only the lines marked "exact min" are proved.  Numerical non-existence = failed
   real-number search with few starts; it says nothing rigorous, and nothing about mod-2^p (2-adic
   border) identities, which may be cheaper.  The N=5 border candidate at u=3 is unverified.
 - Remark (i) of Section 2 (the flattening family is capped at O(M/sqrt u)) and remark (ii) (general
   partitions) are sketches, not proofs.  The HWY/Hurwitz facts are recalled from the literature.

Files: notes.md, thmAprime.py (+_out.txt), exact.py (+_out.txt), search.py, search2.py, search3.py, out_*.txt
